# Input Feature Specification

This document describes the definition, source, computation and preprocessing of the input features that the `03-train_esscls` · `03-train_fitpred` · `04-ablation_esscls` · `04-ablation_fitpred` notebooks feed into the model.
External facts not found in the notebook code (model size, Biopython metric definitions, etc.) were checked against the sources listed at the end of this document.

## Table of Contents

1. [Overview](#1-overview)
2. [Common Conventions](#2-common-conventions)
3. [Variable-Length Inputs](#3-variable-length-inputs) — `esm` · `seq`
4. [Dataset Blocks](#4-dataset-blocks) — `esm_mean` · `prot` · `bio` · `ctx` · `ctxp` · `ctxb`
5. [Final Dataset Standardization](#5-final-dataset-standardization)
6. [Non-Input Columns](#6-non-input-columns)
7. [Inputs Used per Ablation](#7-inputs-used-per-ablation)
8. [Sources](#8-sources)

---

## 1. Overview

| Name | Kind | Shape per sample | Width in concatenated vector | Source file | Processing inside the model |
|---|---|---|---|---|---|
| `seq` | variable-length | `[L]` int8 | 64 | `seq_aa` of the base data | embedding + positional embedding → PMA → `[64]` |
| `esm` | variable-length | `[T, 1152]` fp16 | 1152 | `data/emb_esm_600m.h5` | PMA → `[1152]` |
| `esm_mean` | dataset block | `[1152]` | 1152 | `data/emb_esm_600m_mean.h5` | as is (ablation only; takes the same slot in place of `esm`) |
| `prot` | dataset block | `[1152]` | 1152 | `data/emb_proteome_328m.h5` | as is |
| `bio` | dataset block | `[140]` | 140 | `seq_aa` · `seq_dna` | as is |
| `ctx` | dataset block | `[68]` | 68 | ortholog list (ESM mean space) + reference labels | as is |
| `ctxp` | dataset block | `[68]` | 68 | ortholog list (ProteomeLM space) + reference labels | as is |
| `ctxb` | dataset block | `[68]` | 68 | ortholog list (combined space of both) + reference labels | as is |

- **Dataset** — the fixed-length feature matrix (samples × columns) formed by concatenating the dataset blocks. It resides on the GPU and is rebuilt for every inner fold (because the context features depend on the fold)
- **Concatenated vector** — the vector formed inside the model by concatenating the pooled outputs of `seq` · `esm` and the dataset, in that order. With the full input it has $64 + 1152 + 1152 + 140 + 3 \times 68 = 2712$ dimensions
- **Sample** — one gene. esscls uses 195,354 samples (essential · non-essential), and fitpred uses 177,317 samples (genes with a numeric fitness value)
- Embeddings and sequence features are computed on **all 217,468 proteins, including genes with unknown labels**, and then only the samples in use are selected

---

## 2. Common Conventions

### 2.1 Sample Alignment

Line $i$ of the per-strain base data TSVs (`data/basedata/<index>_<org_id>.tsv`) concatenated in order of their index = row $i$ of each h5 file. The notebooks apply the TSV row numbers `keep` of the samples in use directly to the h5 files as well, and verify by matching `gene_id`.

### 2.2 Per-Strain z-Normalization (`zscore_per_org`)

`esm_mean` · `prot` · `bio` are z-normalized column by column **separately for each strain (proteome)**.

$$\tilde{x}_{ij} = \frac{x_{ij} - \mu_{o(i),j}}{\sigma_{o(i),j} + 10^{-6}}$$

- $x_{ij}$ — the value of column $j$ for gene $i$, $o(i)$ — the strain of gene $i$
- $\mu_{o,j}$ · $\sigma_{o,j}$ — the mean · standard deviation of column $j$ over **all genes** of strain $o$ (including those with unknown labels)
- Meaning: removes genome-level offsets that differ between strains (e.g. genome-wide GC content), leaving only "how the gene compares relative to others within that strain"
- Example: if GC content in strain A has mean 0.65 · standard deviation 0.03, a gene with GC 0.71 → $(0.71 - 0.65)/0.03 = 2.0$

The result is stored as fp16.

---

## 3. Variable-Length Inputs

### 3.1 `esm` — ESM-C 600M Per-Token Embeddings

| Item | Details |
|---|---|
| Model | ESM C 600M (`esmc_600m`, EvolutionaryScale, released 2024-12). A transformer protein language model with 36 layers · width 1152 · 18 attention heads |
| Computation | `02-embed_proteome` Stage A. Each protein sequence is truncated to its first 4,096 residues, fed into ESM-C, and the output `embeddings` are stored per token (bf16 inference → fp16 storage) |
| Shape | `[T, 1152]` per protein, $T$ = number of residues (≤ 4,096) + 2 (including BOS · EOS) |
| Storage | `emb_esm_600m.h5` — `emb` `[total number of tokens, 1152]` · `offsets` `[n + 1]` · `is_special` (1 = BOS/EOS) · `is_truncated` |
| Scale | 71,840,431 tokens in total, about 165.5 GB (fp16). 80 proteins (0.037%) truncated at 4,096 residues |
| Loading | Only the tokens of the samples in use are concatenated in RAM (`esm`, sample $i$ = `esm[esm_off[i]:esm_off[i+1]]`) |

**Processing inside the model — Pooling by Multihead Attention (PMA)**

$$\mathbf{a} = \mathrm{LN}\big(\mathbf{q} + \mathrm{MHA}(\mathbf{q},\ K{=}V{=}\mathbf{E})\big)$$

- $\mathbf{E}$ — token embeddings `[T, 1152]` (padding positions within the batch are excluded by a mask)
- $\mathbf{q}$ — one learnable seed vector `[1152]` (initialized as $0.02 \cdot \mathcal{N}(0, 1)$)
- $\mathrm{MHA}$ — multi-head attention with 8 heads (`N_HEAD`). Since there is only one query, the result is an attention-weighted average of the tokens
- $\mathrm{LN}$ — LayerNorm. The output $\mathbf{a}$ `[1152]` enters the concatenated vector
- Example: initially $\mathbf{q}$ is small, so the attention weights are nearly uniform; $\mathbf{a}$ therefore starts close to the token mean and, during training, concentrates weight on important tokens
- Computation scheme: because the query is a single seed independent of the input, the score of head $h$ is rewritten as $\mathbf{q}_h^\top (W_{k,h}\mathbf{e}_t + \mathbf{b}_{k,h}) = (W_{k,h}^\top \mathbf{q}_h)^\top \mathbf{e}_t + \text{constant}$ and the output as $\sum_t \alpha_t (W_{v,h}\mathbf{e}_t + \mathbf{b}_{v,h}) = W_{v,h}\big(\sum_t \alpha_t \mathbf{e}_t\big) + \big(\sum_t \alpha_t\big)\mathbf{b}_{v,h}$. The constant cancels out in the softmax, so the values are identical to `nn.MultiheadAttention`, and the per-token 1152 → 2304 key · value projection is eliminated

### 3.2 `seq` — Integer-Encoded Sequence

| Item | Details |
|---|---|
| Source | `seq_aa` of the base data (protein sequences from Fitness Browser `aaseqs`) |
| Encoding | 0 = padding, 1 ~ 20 = standard amino acids in the order `ACDEFGHIKLMNPQRSTVWY`, 21 = any character outside the 20 standard ones (X · B · Z · U, etc.) |
| Length | first `MAX_LEN` = 4,096 residues (same cap as the ESM input) |
| Storage | All samples are concatenated into one array `seq_flat`, with boundaries given by `seq_off`. During training the whole array is placed on the GPU, and each batch `[B, L]` is built on the GPU by index slicing (padding 0) |

**Processing inside the model** — `SeqBranch`

$$\mathbf{s} = \mathrm{PMA}\big(\mathrm{Emb}(x_{1:L}) + \mathbf{P}_{1:L}\big)$$

- $\mathrm{Emb}$ — embedding of the 22-token vocabulary → 64 dimensions (`SEQ_DIM`) (index 0 = fixed padding vector)
- $\mathbf{P}$ — learnable positional embedding `[4096, 64]`
- $\mathrm{PMA}$ — the same pooling as in Section 3.1, with 4 heads (`SEQ_HEAD`). Output $\mathbf{s}$ `[64]`
- Example: `MKT` → `[11, 9, 17]` → (amino acid vector + position vector) at each position → pooling → one 64-dimensional vector

---

## 4. Dataset Blocks

### 4.1 `esm_mean` — ESM-C Token-Mean Embedding (1152 columns, ablation only)

| Item | Details |
|---|---|
| Computation | Per-protein arithmetic mean of the token embeddings from Section 3.1 (including BOS · EOS) |
| Storage | `emb_esm_600m_mean.h5` `[217468, 1152]` |
| Preprocessing | per-strain z-normalization → selection of samples in use → fp16 |
| Usage | Used in place of `esm` (tokens + PMA) in the `poolesm` · `poolesm-noseq` · `poolesm-noder` · `poolesm-prot` ablations. In the 03 notebooks it is not an input and is used only as the similarity space of `ctx` · `ctxb` |

### 4.2 `prot` — ProteomeLM-L Contextual Embedding (1152 columns)

| Item | Details |
|---|---|
| Model | ProteomeLM-L (`Bitbol-Lab/ProteomeLM-L`, 328M parameters). A transformer that takes an entire proteome as input at once and is trained to reconstruct masked protein embeddings from the proteome context |
| Input | `[n_org, 1152]`, the per-protein ESM-C 600M token-mean vectors (including BOS · EOS, the same vectors as in Section 4.1) laid out for **all proteins of a single strain** |
| Output | `hidden_states[8]` — the output of the 8th of 18 blocks (`hidden_states[0]` = input embedding). `[1152]` per protein |
| Computation | `02-embed_proteome` Stage B. The 46 strains are passed through one at a time |
| Preprocessing | per-strain z-normalization (Section 2.2) → selection of samples in use → fp16 |
| Meaning | A protein representation reflecting its relationships with the other proteins in the same genome |

The same raw vectors (before normalization) are used in two more places.

- **Ortholog search** — the similarity space of `ctxp` · `ctxb` (Section 4.4)
- **Strain proteome representation** `ORG_VEC` — the per-strain mean of the raw `prot`. Used only for selecting the validation genus; it is not a model input

### 4.3 `bio` — Biopython Sequence Features (140 columns)

72 protein-sequence feature columns + 68 CDS (protein-coding sequence) feature columns. Computed on all 217,468 genes → per-strain z-normalization → selection of samples in use.

#### (1) 12 Global Physicochemical Feature Columns — `Bio.SeqUtils.ProtParam.ProteinAnalysis`

Before computation, characters other than the 20 standard ones are removed. The remaining sequence is denoted $s$ and the original length $L$.

| # | Feature | Computation | Scaling |
|---|---|---|---|
| 1 | Length | $\log(1 + L)$ — $L$ is the original length including non-standard characters | — |
| 2 | Molecular weight | `molecular_weight()` (Da) | ÷ 10⁴ |
| 3 | Isoelectric point | `isoelectric_point()` — Biopython `IsoelectricPoint` module | — |
| 4 | Aromaticity | `aromaticity()` — Lobry (1994), relative frequency of Phe + Trp + Tyr | — |
| 5 | Instability index | `instability_index()` — Guruprasad et al. (1990). A value above 40 is regarded as unstable (short half-life) | ÷ 100 |
| 6 | GRAVY | `gravy()` — mean of the Kyte · Doolittle (1982) hydropathy scale (Grand Average of Hydropathy) | — |
| 7 | Helix-propensity fraction | `secondary_structure_fraction()[0]` — fraction of E · M · A · L · K | — |
| 8 | Turn-propensity fraction | `secondary_structure_fraction()[1]` — fraction of N · P · G · S · D | — |
| 9 | Sheet-propensity fraction | `secondary_structure_fraction()[2]` — fraction of V · I · Y · F · W · L · T | — |
| 10 | Charge | `charge_at_pH(7.0)` | ÷ 10 |
| 11 | Molar extinction coefficient (reduced) | `molar_extinction_coefficient()[0]` — assumes all cysteines are in the reduced form | ÷ 10⁴ |
| 12 | Molar extinction coefficient (oxidized) | `molar_extinction_coefficient()[1]` — assumes cysteines form disulfide bonds (cystine) | ÷ 10⁴ |

- 7 ~ 9 are **structural propensities inferred from amino acid composition, not actual secondary structure**. L is counted in both helix and sheet, so the three values do not sum to 1
- **Version caveat** — before Biopython 1.82, this method returned the values in the order (sheet, turn, helix) (a fact stated in the 1.82 documentation). The code assumes the order (helix, turn, sheet) via `h, t, e = ...`. **Assumption**: the model training environment was taken to be 1.82 or later. This is because the installation log of the 02 notebook records biopython 1.88; the version in the execution environment of the 03 · 04 notebooks was not checked directly
- If there are fewer than 5 standard residues, all 72 columns are set to 0 and only column 1 (length) is filled

#### (2) 60 Amino Acid Composition Columns

| Group | Columns | Computation (for each amino acid $a$) |
|---|---|---|
| Overall composition | 20 | $\mathrm{count}_s(a) / \lvert s \rvert$ |
| N-terminal composition | 20 | $\mathrm{count}_{s[:20]}(a) / \min(20, \lvert s \rvert)$ — first 20 residues |
| C-terminal composition | 20 | $\mathrm{count}_{s[-20:]}(a) / \min(20, \lvert s \rvert)$ — last 20 residues |

- Example: if the first 20 residues of $s$ contain 1 M and 3 K, the N-terminal composition has M = 0.05, K = 0.15

#### (3) 68 CDS Feature Columns — `seq_dna`

`seq_dna` is the gene sequence cut out of the genome (reverse-complemented for the reverse strand, stop codon included).

| Group | Columns | Computation |
|---|---|---|
| Codon frequency | 64 | Codons in the reading frame are counted 3 bases at a time from the start and converted into relative frequencies of the 64 codons in `TCAG` order. Codons containing characters other than A · C · G · T are excluded. Stop codons are also counted |
| GC | 1 | G + C fraction of the whole sequence |
| GC1 · GC2 · GC3 | 3 | G + C fraction at codon positions 1 · 2 · 3 |

- If there are fewer than 3 codons, all 68 columns are set to 0
- Example: `ATGGCCTAA` → codons `ATG` · `GCC` · `TAA` at 1/3 each, GC = 4/9, GC1 = 1/3 (G among A · G · T), GC2 = 1/3 (C among T · C · A), GC3 = 2/3 (G · C among G · C · A)

### 4.4 `ctx` · `ctxp` · `ctxb` — Cross-Strain Context Features (68 columns each)

These summarize the **labels of the most similar proteins (ortholog candidates)** in other strains and the **extent to which similar proteins exist** within the gene's own genome. The three blocks share the same formulas and **differ only in the similarity space**.

| Block | Similarity space (`CTX_SPACE`) | Embedding |
|---|---|---|
| `ctx` | `esm` | ESM-C token mean (raw, Section 4.1) |
| `ctxp` | `prot` | ProteomeLM (raw, Section 4.2) |
| `ctxb` | `both` | the two whitened embeddings concatenated and divided by $\sqrt{2}$ |

#### (1) Ortholog List `ORT` (fold-independent, computed once)

1. **Whitening** — the embeddings are standardized per dimension and then L2-normalized. The inner product of two vectors then equals their cosine similarity
2. **Best hit per strain** — for each of all 217,468 proteins, the index and cosine of the protein with the highest cosine in each of the 46 strains are stored (`sim`, `idx` `[217468, 46]`). The protein itself is excluded
3. **Top cosines within own genome** — the top 20 cosines within the same strain are stored in descending order (`para` `[217468, 20]`, excluding the protein itself)

#### (2) Reference Strains and Top 20 (computed per fold)

For a query sample $i$, a **reference strain** is a strain that belongs to none of the following three.

- test genus (the test genus of the outer fold)
- validation genus (the validation genus of the inner fold)
- the genus of sample $i$ itself

Only the best hits from reference strains are sorted by descending cosine, and the top 20 (`max(TOPK)`) are used. The reference protein at rank $r = 1, \dots, 20$ is denoted $J_r$ and its cosine $S_r$. If there are fewer than 20 reference strains, the remaining slots are left invalid.

Because of this rule, the labels of the test · validation genera never enter the input of any sample in that fold. In the final model (no test · validation genus), only the sample's own genus is excluded.

#### (3) Reference Labels

| Symbol | Definition |
|---|---|
| $f_r$ | `fit_mean` of reference protein $J_r$ ($-\infty$ if essential, NaN if unknown) |
| $\mathrm{ess}_r$ | $J_r$ is essential ($f_r = -\infty$) |
| $\mathrm{fin}_r$ | $J_r$ has a numeric fitness value |
| $\mathrm{strong}_r$ | $\mathrm{fin}_r$ and $\lvert \mathrm{comb\_t} \rvert \ge 2$ (`T_EVAL`) |
| $q_r$ | **within-strain fitness rank** of $J_r$, in $[-1, 1]$ |
| $\mathrm{rbh}_r$ | reciprocal best hit: the best hit of $i$ in the strain of $J_r$ is $J_r$, and at the same time the best hit of $J_r$ in the strain of $i$ is $i$ |

The within-strain rank is computed separately for each strain, using only genes with numeric fitness values.

$$q = \frac{\mathrm{rank} + 0.5}{n_{\mathrm{org}}} \times 2 - 1$$

- $\mathrm{rank}$ — ascending rank within the strain (starting from 0), $n_{\mathrm{org}}$ — the number of genes with numeric fitness values in that strain
- Example: out of 4,000 genes, the lowest fitness → $q \approx -1.0$, the median → $q \approx 0$, the highest fitness → $q \approx +1.0$
- Meaning: removes fitness scales that differ between strains (experimental conditions · noise levels), leaving only the relative position

The similarity weight is $w_r = \exp\big(5\,(S_r - 1)\big)$. Example: $S = 1.0 \to 1.00$, $S = 0.8 \to 0.37$, $S = 0.5 \to 0.08$ — the more similar the reference, the more strongly it is weighted.

#### (4) Composition of the 68 Columns

Below, "top $k$" means the valid references at ranks 1 ~ $k$, with $k \in$ `TOPK` = {1, 3, 5, 10, 20}. If there is nothing to average or take a ratio over, the value is 0, and any NaN · ±∞ in the computed results is also replaced with 0.

**`cons_*` — degree of conservation (9 columns, no labels used)**

| Column | Computation |
|---|---|
| `cons_simtop{k}` (5) | mean of the top $k$ cosines |
| `cons_gap1_3` | $S_1 - S_3$ — how much rank 1 stands out from rank 3 |
| `cons_gap1_10` | $S_1 - S_{10}$ |
| `cons_sim_sd` | standard deviation of the top 20 cosines (invalid slots are 0) |
| `cons_nref` | number of valid references / 20 |

**`otf_*` — reference fitness summary, original scale (25 columns)**

| Column | Computation |
|---|---|
| `otf_mean{k}` (5) | mean of $f_r$ over the references in the top $k$ that are $\mathrm{fin}$ |
| `otf_essfrac{k}` (5) | fraction essential among the references in the top $k$ that have a label ($\mathrm{fin}$ or $\mathrm{ess}$) |
| `otf_wmean{k}` (5) | $w_r$-weighted mean $\sum w_r f_r / \sum w_r$ over the references in the top $k$ that are $\mathrm{fin}$ |
| `otf_mean_strong10` | mean of $f_r$ over the references in the top 10 that are $\mathrm{strong}$ |
| `otf_nstrong10` | number of references in the top 10 that are $\mathrm{strong}$ / 10 |
| `otf_min10` · `otf_max10` | minimum · maximum of $f_r$ over the references in the top 10 that are $\mathrm{fin}$ |
| `otf_sd10` | standard deviation of $f_r$ over the top 10 (slots that are not $\mathrm{fin}$ are 0) |
| `otf_best_fit` | $f_1$ of the rank-1 reference (when $\mathrm{fin}$) |
| `otf_best_sim` | cosine $S_1$ of the rank-1 reference |
| `otf_best_isess` | 1 if the rank-1 reference is essential |
| `otf_frac_lt_m1` · `otf_frac_gt_p1` | number of references in the top 10 with $f_r < -1$ · $f_r > 1$ / 10 |

**`otr_*` — reference fitness summary, within-strain rank scale (12 columns)**

| Column | Computation |
|---|---|
| `otr_mean{k}` (5) | mean of $q_r$ over the top $k$ |
| `otr_wmean{k}` (5) | $w_r$-weighted mean of $q_r$ over the top $k$ |
| `otr_min10` | minimum of $q_r$ over the top 10 |
| `otr_frac_low20` | number of references in the top 20 with $q_r < -0.5$ (bottom 25% of the strain) / 20 |

**`rbh_*` — reciprocal best hits only (14 columns)**

| Column | Computation |
|---|---|
| `rbh_frac` | number of reciprocal best hits in the top 20 / 20 |
| `rbh_sim_mean` | mean cosine of the reciprocal-best-hit references |
| `rbh_n{k}` ($k$ = 5, 10, 20) | number of reciprocal best hits with fitness values in the top $k$ / $k$ |
| `rbh_fit{k}` | mean of $f_r$ over those references |
| `rbh_rank{k}` | mean of $q_r$ over those references |
| `rbh_essfrac{k}` | fraction essential among the labeled references among the reciprocal best hits in the top $k$ |

**`para_*` — similar proteins within own genome (8 columns, no labels used)**

Computed from the top 20 within-genome cosines $P_1 \ge \dots \ge P_{20}$ from step (1)-3.

| Column | Computation |
|---|---|
| `para_simtop{k}` (5) | mean of $P_1, \dots, P_k$ |
| `para_gap1_3` · `para_gap1_10` | $P_1 - P_3$ · $P_1 - P_{10}$ |
| `para_sim_sd` | standard deviation of $P_1, \dots, P_{20}$ |

**Assumption**: the design intent is presumed to rest on the hypothesis that "if there are many similar proteins (paralog candidates), the function may be compensated even when the gene is lost" (an inference; the code does not state a rationale). The **labels** of paralogs are not used, because if the query is a test-genus sample, its paralogs' labels also belong to the test genus.

**Example** — if the top 3 references of the query gene are (cosine 0.92, essential) · (0.85, $f = -2.1$) · (0.60, $f = 0.3$):
`cons_simtop3` = 0.79, `otf_essfrac3` = 1/3, `otf_mean3` = $(-2.1 + 0.3)/2 = -0.9$, `otf_best_isess` = 1, `otf_best_sim` = 0.92.

---

## 5. Final Dataset Standardization

`build_dataset` concatenates the required blocks and then standardizes all samples using the column-wise mean · standard deviation of **the training samples of that inner fold**.

$$z_{ij} = \frac{x_{ij} - \bar{x}_j^{\,\mathrm{train}}}{s_j^{\,\mathrm{train}} + 10^{-6}}$$

- $\bar{x}_j^{\,\mathrm{train}}$ · $s_j^{\,\mathrm{train}}$ — mean · standard deviation of column $j$ over the training samples (validation · test samples are not used)
- The final model uses all samples as the reference
- NaN is replaced with 0, and the result is loaded as an fp16 GPU tensor

Per-strain z-normalization (Section 2.2) uses label-free within-strain statistics, whereas this standardization is a global rescaling redone for every fold.

---

## 6. Non-Input Columns

| Column | Usage |
|---|---|
| `fit_mean` | Source of the labels (esscls: $-\infty$ → essential, fitpred: $> 0$ → high). The values of reference proteins are used in the context features |
| `comb_t` | fitpred training-sample weight clip(\|comb_t\|) and the criterion for early-stopping validation samples (\|comb_t\| ≥ 2). The values of reference proteins are used in `otf_*strong*` |
| `org_id` · `genus` | fold splitting · sample weights · per-strain metrics |
| `ORG_VEC` | strain proteome representation for validation genus selection |

A sample's own `fit_mean` · `comb_t` never enters any input feature.

---

## 7. Inputs Used per Ablation

| ablation | `seq` | `esm` | `esm_mean` | `prot` | `bio` | `ctx` · `ctxp` · `ctxb` | Concatenated vector width |
|---|:-:|:-:|:-:|:-:|:-:|:-:|--:|
| `full` (03 notebooks) | ✓ | ✓ | | ✓ | ✓ | ✓ | 2712 |
| `noctx` | ✓ | ✓ | | ✓ | ✓ | | 2508 |
| `nobio` | ✓ | ✓ | | ✓ | | ✓ | 2572 |
| `noesm` | ✓ | | | ✓ | ✓ | ✓ | 1560 |
| `noseq` | | ✓ | | ✓ | ✓ | ✓ | 2648 |
| `noder` | ✓ | ✓ | | ✓ | | | 2368 |
| `noemb` | ✓ | | | | ✓ | ✓ | 408 |
| `seq-prot` | ✓ | | | ✓ | | | 1216 |
| `esm-prot` | | ✓ | | ✓ | | | 2304 |
| `prot` | | | | ✓ | | | 1152 |
| `poolesm` | ✓ | | ✓ | ✓ | ✓ | ✓ | 2712 |
| `poolesm-noseq` | | | ✓ | ✓ | ✓ | ✓ | 2648 |
| `poolesm-noder` | ✓ | | ✓ | ✓ | | | 2368 |
| `poolesm-prot` | | | ✓ | ✓ | | | 2304 |

- The row order matches the execution order of `ABLATIONS` in the 04 notebooks (`full` reads the 03 notebook results as the baseline for comparison)
- The "der" in `noder` · `poolesm-noder` stands for derived features, i.e. `bio` · `ctx` × 3 computed from the raw sequences · embeddings
- The concatenation order inside the model also follows the column order of the table — `seq` → `esm` → dataset blocks (`esm_mean` → `prot` → `bio` → `ctx` → `ctxp` → `ctxb`), including only those present in that ablation

---

## 8. Sources

External references

- ESM Cambrian (model introduction), EvolutionaryScale blog, 2024-12-04 — https://www.evolutionaryscale.ai/blog/esm-cambrian
- EvolutionaryScale/esmc-600m-2024-12 model card, Hugging Face — https://huggingface.co/EvolutionaryScale/esmc-600m-2024-12
- Malbranke C., Zalaffi G. P., Bitbol A.-F., *ProteomeLM: A proteome-scale language model enables accurate and rapid prediction of protein–protein interactions and gene essentiality across taxa*, bioRxiv 2025 (doi 10.1101/2025.08.01.668221) — https://www.biorxiv.org/content/10.1101/2025.08.01.668221v2 · PNAS published version — https://www.pnas.org/doi/abs/10.1073/pnas.2524201123
- Bitbol-Lab/ProteomeLM official repository (model size table: ProteomeLM-L 328M), GitHub — https://github.com/Bitbol-Lab/ProteomeLM
- Bio.SeqUtils.ProtParam module documentation (definitions of each metric · cited references), Biopython — https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html
- The same documentation, version 1.82 (record of the fix to the return order of secondary structure fractions), Biopython — https://biopython.org/docs/1.82/api/Bio.SeqUtils.ProtParam.html
- Lee J. et al., *Set Transformer: A Framework for Attention-based Permutation-Invariant Neural Networks*, ICML 2019, PMLR 97:3744–3753 — https://proceedings.mlr.press/v97/lee19d.html
- Price M. N. et al., *Mutant phenotypes for thousands of bacterial genes of unknown function*, Nature 557:503–509, 2018 (Fitness Browser) — https://www.nature.com/articles/s41586-018-0124-0

Project code

- `01-build_basedata.ipynb` — definitions of `fit_mean` · `comb_t` · `seq_aa` · `seq_dna`
- `02-embed_proteome.ipynb` — computation of ESM-C token embeddings · token means · ProteomeLM `hidden_states[8]`; token counts · truncation counts · number of hidden layers
- `03-train_esscls.ipynb` · `03-train_fitpred.ipynb` · `04-ablation_*.ipynb` — execution of each stage
- `proteomefit/features.py` · `proteomefit/models.py` — feature formulas · preprocessing · processing inside the model
