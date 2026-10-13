# ProteomeFit

Experiment code for predicting the **essentiality** and **knockout fitness** of bacterial genes from protein language model embeddings and **proteome context**.
Generalization to unseen genera is evaluated with genus-level nested cross-validation.

- Data: Fitness Browser (RB-TnSeq), 46 bacterial strains · 28 genera · 217,468 protein-coding genes
- Task 1 — **essential gene classification** (`esscls`): 195,354 samples, 18,037 essential (9.23%)
- Task 2 — **fitness regression** (`fitpred`): 177,317 samples with a numeric fitness value; the target is `fit_mean` z-normalized within each strain
- Inputs: ESM C 600M per-token embeddings · ProteomeLM-L proteome-context embeddings · Biopython sequence features · cross-strain context features (summaries of ortholog labels)

> [!TIP]
> **Try it online: [https://proteomefit.duckdns.org](https://proteomefit.duckdns.org)**
>
> Upload the genes of one bacterial strain (GenBank, CDS FASTA or complete-genome FASTA) and get, for every gene, an essentiality call and, for non-essential genes, the predicted knockout fitness (within-strain z-score) with its uncertainty, as a downloadable CSV. The service runs the final models from `results/esscls_ablation/nobio` (classification) and `results/fitpred_ablation/noseq` (regression). Uploaded files and results are not stored on the server.

## Main results

Values are the mean ± standard deviation across genera of the per-genus performance (each genus = mean of its 3 inner folds), taken from the `mean` and `std` rows of `results/*/metrics_total.csv`.

**Essential gene classification** — 26 genera (2 genera without essential genes are excluded)

| Result folder | Dropout | MCC | ROC AUC | PR AUC |
|---|---|---|---|---|
| `esscls_10` | 0.1 | 0.705 ± 0.077 | 0.932 ± 0.042 | 0.798 ± 0.075 |
| `esscls_30` | 0.3 | 0.712 ± 0.078 | 0.936 ± 0.040 | 0.811 ± 0.074 |
| `esscls_50` | 0.5 | 0.708 ± 0.076 | 0.938 ± 0.040 | 0.820 ± 0.073 |

**Fitness regression** — 28 genera. `_sig` = test-genus genes with |comb_t| ≥ 2, whose measured fitness is distinguishable from 0 (main metric; the same genes for every model); `_all` = all test genes

| Result folder | Dropout | Spearman_sig | CCC_sig | Spearman_all |
|---|---|---|---|---|
| `fitpred_10` | 0.1 | 0.381 ± 0.130 | 0.328 ± 0.118 | 0.118 ± 0.056 |
| `fitpred_30` | 0.3 | 0.380 ± 0.129 | 0.322 ± 0.119 | 0.118 ± 0.055 |
| `fitpred_50` | 0.5 | 0.375 ± 0.127 | 0.315 ± 0.116 | 0.116 ± 0.053 |

**Input ablation** — difference in the genus-mean metric relative to the full-input model (classification `esscls_30`, regression `fitpred_10`)

| Ablation | Inputs | ΔMCC (classification) | ΔSpearman_sig (regression) |
|---|---|---|---|
| `noctx` | seq + esm + prot + bio | −0.025 | −0.079 |
| `nobio` | seq + esm + prot + ctx × 3 | +0.002 | −0.002 |
| `noesm` | seq + prot + bio + ctx × 3 | −0.003 | −0.004 |
| `noseq` | esm + prot + bio + ctx × 3 | −0.001 | +0.002 |
| `noder` | seq + esm + prot | −0.028 | −0.079 |
| `noemb` | seq + bio + ctx × 3 | −0.054 | −0.009 |
| `seq-prot` | seq + prot | −0.049 | −0.116 |
| `esm-prot` | esm + prot | −0.026 | −0.071 |
| `prot` | prot | −0.046 | −0.079 |
| `poolesm` | seq + esm_mean + prot + bio + ctx × 3 | −0.009 | −0.003 |
| `poolesm-noseq` | esm_mean + prot + bio + ctx × 3 | −0.013 | +0.001 |
| `poolesm-noder` | seq + esm_mean + prot | −0.039 | −0.096 |
| `poolesm-prot` | esm_mean + prot | −0.035 | −0.081 |

seq = integer-encoded sequence, esm = ESM per-token embeddings, esm_mean = mean of ESM token embeddings, prot = ProteomeLM embedding, bio = Biopython features, ctx × 3 = the three context-feature blocks.
Full tables: `results/esscls_ablation/ablation_summary.csv` · `results/fitpred_ablation/ablation_summary.csv`.

## Repository layout

```
proteomefit/
├── README.md
├── data/
│   ├── basedata/                 # 46 per-strain base data TSVs (01_acidovorax_3H11.tsv … 46_WCS417.tsv)
│   └── refer_list.md             # publications and supplementary tables used by notebook 01
├── results/                      # experiment outputs (final models model_final.pt included; per-fold model.pt excluded)
│   ├── esscls_{10,30,50}/        # essential gene classification, suffix = dropout x 100
│   ├── fitpred_{10,30,50}/       # fitness regression
│   ├── esscls_ablation/          # 13 classification ablations
│   ├── fitpred_ablation/         # 13 regression ablations
│   └── calib_conf/               # fitness-confidence calibration (conf_sign · conf_fit)
└── src/
    ├── 01-build_basedata.ipynb   # Fitness Browser -> per-strain base data TSVs
    ├── 02-embed_proteome.ipynb   # ESM C 600M · ProteomeLM-L embedding h5 files
    ├── 03-train_esscls.ipynb     # essential gene classification (nested CV + final model)
    ├── 03-train_fitpred.ipynb    # fitness regression
    ├── 04-ablation_esscls.ipynb  # classification input ablation
    ├── 04-ablation_fitpred.ipynb # regression input ablation
    ├── 05-calib_conf.ipynb       # calibration of fitness confidences
    ├── input_features.md         # input feature specification (definitions · formulas · preprocessing)
    └── proteomefit/              # shared modules imported by the notebooks
```

## Pipeline

| Step | Notebook | Output |
|---|---|---|
| 1 | `01-build_basedata` | `data/basedata/<index>_<org_id>.tsv` — per strain: coding genes · sequences · summary fitness values in no-stress conditions · combined statistics · essentiality |
| 2 | `02-embed_proteome` | `data/emb_esm_600m.h5` (per token) · `data/emb_esm_600m_mean.h5` (token mean) · `data/emb_proteome_328m.h5` (ProteomeLM `hidden_states[8]`) |
| 3 | `03-train_esscls` · `03-train_fitpred` | `results/<task>_<dropout>/` — per-fold models · metrics · figures, `metrics_total.csv`, final model |
| 4 | `04-ablation_*` | `results/<task>_ablation/<ablation>/` · `ablation_summary.csv` |
| 5 | `05-calib_conf` | `results/calib_conf/calibration.json` — calibrations of the regression mean · SD into `conf_sign` (probability that the measured fitness lies on the predicted side of the strain average) and `conf_fit` (probability that it lies within mean ± SD) — gradient-boosting probability tables fitted on held-out-genus predictions of all genes |

### Base data (`data/basedata/`)

There is one TSV per strain; the leading index in each file name is the strain order in notebook 01. **Concatenating the files in index order gives the row order of the embedding h5 files**, and `proteomefit.data.load_basedata` reads them in that order.

| Column | Content |
|---|---|
| `org_id` · `org_name` | strain ID · scientific name |
| `gene_id` · `gene_name` | locus ID · gene name |
| `seq_dna` · `seq_aa` | gene DNA sequence (reverse complement for minus-strand genes, stop codon included) · protein sequence |
| `fit_mean` · `fit_median` | mean · median fitness over no-stress experiments (essential = `-inf`, unknown = empty) |
| `comb_se` · `comb_t` | standard error · t-like statistic combining replicate experiments |
| `list_fit` · `list_se` · `list_t` | per-experiment raw values |

### Supplementary essentiality tables

For strains not covered by the Fitness Browser file `essential_proteins.tab`, notebook 01 reads essential gene lists from the supplementary tables of the original publications (Liu21 · Trotter21 · Helmann19 · Helmann21 · Georgoulis21). These files are not included for copyright reasons; the publications and file names are listed in [`data/refer_list.md`](data/refer_list.md). Put the downloaded files in `data/refer/` to rerun notebook 01.

### Files not included in the repository

| File | Size | How to obtain |
|---|---|---|
| `data/refer/` (supplementary tables · papers) | ~52 MB | supplementary material of the publications in `data/refer_list.md` |
| `data/fitness_browser/` (`feba.db` · `aaseqs` · `essential_proteins.tab`) | ~9.7 GB | [Fitness Browser](https://fit.genomics.lbl.gov/) data (2024-02 snapshot) |
| `data/emb_esm_600m.h5` | ~165.6 GB | notebook 02, Stage A |
| `data/emb_esm_600m_mean.h5` · `data/emb_proteome_328m.h5` | ~0.52 GB each | notebook 02, Stage B |
| `results/**/fold*/inner*/model.pt` (per-fold models) | ~63.0 GB in total | notebooks 03 · 04 |

## Method summary

- **Model** — the integer-encoded sequence and the ESM token sequence are each reduced to one vector by Pooling by Multihead Attention (PMA) and concatenated with a fixed-length dataset (ProteomeLM · Biopython · context features), 2,712 dimensions for the full input.
  - Classification: multi-layer perceptron (MLP) `2712 → 640 → 320 → 1`; loss = binary cross-entropy with positive-class weighting and inverse taxon-size sample weights
  - Regression: 1,000 mini-MLPs, each using randomly chosen input blocks and columns (with a per-branch bootstrap); loss = per-branch weighted 1 − concordance correlation coefficient (CCC). The mean of the branch predictions is the prediction and their standard deviation is the uncertainty
- **Context features** — for each query gene, the most similar protein (ortholog candidate) is found in every strain, and only the labels from strains outside the test, validation and query's own genus are summarized over the top k hits. Labels of the test genus therefore never enter any input of that fold
- **Evaluation** — outer folds = leave-one-genus-out (28), inner folds = 3 validation genera whose extrapolation distance is similar to that of the test genus. Metrics are computed per strain and averaged. For regression, the main metric uses only genes with |comb_t| ≥ 2, because the measured fitness of the other genes is not distinguishable from 0

See [`src/input_features.md`](src/input_features.md) for exact feature definitions.

## Shared modules (`src/proteomefit/`)

| Module | Content |
|---|---|
| `config.py` | hyperparameters (`Config`) · global torch settings |
| `data.py` | loading per-strain TSVs · embedding h5 files · ESM tokens (GPU-resident) |
| `features.py` | Biopython features · integer sequences · ortholog lists · context features · dataset |
| `folds.py` | genus-level nested cross-validation splits |
| `models.py` | `PMA` · `SeqBranch` · `EssNet` (classification) · `FitNet` (regression) |
| `batching.py` | length-aware batching · GPU-resident batch feeding · inference |
| `metrics.py` · `plots.py` | metrics · figures |
| `train.py` | training loop · losses · sample weights |
| `tasks.py` | task definitions `EssTask` · `FitTask` (labels · loss · early-stopping metric · evaluation outputs) |
| `cv.py` | nested cross-validation runner · final model · result summaries |

## Running

Run the notebooks with `src/` as the working directory (paths are relative to `../data` and `../results`).

- Notebooks 03 · 04 keep all ESM token embeddings (bf16, ~131–145 GiB) in GPU memory and therefore need a GPU with ~180 GB of VRAM (e.g. NVIDIA B200)
- Main packages: `torch` · `esm` · `proteomelm` · `transformers` · `h5py` · `biopython` · `scikit-learn` · `scipy` · `pandas` · `matplotlib` · `tqdm` · `openpyxl`
- With `resume=True` (the `Config` default), folds and final models whose outputs already exist are not retrained

## References

- Price M. N. et al., *Mutant phenotypes for thousands of bacterial genes of unknown function*, Nature 557:503–509, 2018 — https://www.nature.com/articles/s41586-018-0124-0
- Fitness Browser — https://fit.genomics.lbl.gov/
- EvolutionaryScale, *ESM Cambrian*, 2024-12 — https://www.evolutionaryscale.ai/blog/esm-cambrian
- Malbranke C., Zalaffi G. P., Bitbol A.-F., *ProteomeLM: A proteome-scale language model allowing fast prediction of protein-protein interactions and gene essentiality across taxa*, bioRxiv 2025 — https://www.biorxiv.org/content/10.1101/2025.08.01.668221 · code: https://github.com/Bitbol-Lab/ProteomeLM
- Lee J. et al., *Set Transformer*, ICML 2019 — https://proceedings.mlr.press/v97/lee19d.html
- Biopython `Bio.SeqUtils.ProtParam` — https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html

## Author

Seongbo Heo — Korea Biofoundry, Korea Research Institute of Bioscience and Biotechnology (KRIBB) (hsb@kribb.re.kr)

## Acknowledgments

This research was funded by the Bio & Medical Technology Development Program (grant number 2018M3A9H3024746 and RS-2024-00335300) of the National Research Foundation, funded by the Ministry of Science and ICT of the Republic of Korea. This research was also funded by the Korea Research Institute of Bioscience and Biotechnology (KRIBB) Research Initiative Program (KGM5402423)
