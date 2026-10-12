# FB feba.db — Publications (pubId), recovery media, and essentiality

This document compiles every pubId in feba.db (2024-02 snapshot, 46 bacteria + 2 archaea = 48 organisms), the linked strains, the evidence for each strain's recovery (library-construction) medium, and the sources and counts of essentiality calls. **Recovery media and essentiality are recorded only from the original sources (paper Methods/supplementary tables) or from values observed directly in feba.db, not by inference** (items that could not be verified are explicitly marked). The scope is the 46 bacterial strains; the 2 archaeal strains (Methanococcus) are excluded.

> **File naming rule (refer/)**: `{pubId}_paper.pdf` (main text) · `Price18_SI.pdf`·`Price18_suppl.xlsx` (Price18 supplementary material) · `{pubId}_{orgId}_essential.{csv,xlsx}` (per-strain essentiality) · `{pubId}_ST#_*.xlsx` (multi-strain supplementary tables). This index file is kept separately at `data/refer_list.md`. The files in `refer/` are not included in this repository; obtain them from the publications listed below and place them in `data/refer/`.

---

## Table 1. Publications (pubId), links, and linked strains

### Table 1A. All pubIds in feba.db (20)

- `URL(feba.db)` = the link in the feba.db `Publication` table (many are bioRxiv preprints). `Published version (verified)` = DOI of the confirmed final published version (only where it differs from feba.db).
- Bold pubId = the **library-construction (recovery-medium evidence) or essentiality-call paper** for that strain. The rest are condition-specific studies that reused that library.

| pubId | URL (feba.db) | Published version (verified) | Linked strains (orgId) |
|---|---|---|---|
| **Wetmore15** | https://doi.org/10.1128/mBio.00306-15 | (same) | Keio, MR1, Phaeo, psRCH2, SB2B |
| **Melnyk15** | https://doi.org/10.1128/mBio.00233-15 | (same) | PS |
| Price17 | https://doi.org/10.1371/journal.pgen.1007147 | (same) | BFirm, HerbieS, Marino, Miya, Phaeo, Smeli, psRCH2, pseudo6_N2E2 |
| Rand17 | https://doi.org/10.1038/s41564-017-0028-z | (same) | Putida |
| Thorgersen17 | https://doi.org/10.3389/fmicb.2017.01529 | (same) | psRCH2 |
| **Price18** | http://genomics.lbl.gov/supplemental/bigfit/ | Nature 10.1038/s41586-018-0124-0 | **32 strains** (all Price18-series strains in Table 2 below) |
| Thompson18 | https://doi.org/10.1101/450254 (bioRxiv) | mBio 10.1128/mbio.02577-18 | Putida |
| **Price18b** | http://dx.doi.org/10.1128/msystems.00297-18 | (same) | BFirm, Burk376, Koxy, Putida, WCS417 |
| **Liu21** | https://linkinghub.elsevier.com/retrieve/pii/S2211-1247(21)00103-0 | Cell Rep 10.1016/j.celrep.2021.108789 | **Btheta**, Cola, Pedo557 |
| **Helmann19** | http://dx.doi.org/10.1073/pnas.1908858116 | (same) | **SyringaeB728a** |
| **Helmann19b** | https://doi.org/10.1101/684605 (bioRxiv) | (published version not confirmed) | SyringaeB728a, SyringaeB728a_mexBdelta |
| Incha20 | https://doi.org/10.1016/j.mec.2019.e00119 | (same) | Putida |
| Thompson20 | https://www.biorxiv.org/content/10.1101/2020.07.04.188060v1 | (published version not confirmed) | Putida |
| Garber20 | https://doi.org/10.1101/2020.09.30.321588 (bioRxiv) | (published version not confirmed) | Putida |
| **Helmann21** | https://doi.org/10.1101/2021.09.16.460530 (bioRxiv) | Front Microbiol 10.3389/fmicb.2022.778927 (PMC8821946) | **Dda3937, Ddia6719, DdiaME23** |
| **McCausland21** | https://doi.org/10.1101/2021.08.13.456315 (bioRxiv) | mSystems 10.1128/msystems.01037-21 (PMC8788322) | **Magneto** |
| **Trotter21** | https://www.biorxiv.org/content/10.1101/2021.01.13.426591v1 | Front Microbiol 10.3389/fmicb.2023.1095191 (PMC10102598) | **DvH** |
| Price21a | https://doi.org/10.1371/journal.pgen.1009342 | (same) | Korea |
| Price21b | https://doi.org/10.1101/2021.11.02.466981 (bioRxiv) | PLoS Genet 10.1371/journal.pgen.1010156 (PMC9007349, 2022) — title changed on publication (preprint "GapMind for Carbon Sources: Automated annotations of catabolic pathways" → published version "Filling gaps in bacterial catabolic pathways with computation and high-throughput genetics"). | ANA3, Burk376, Cup4G11, Dino, HerbieS, Phaeo, Putida, WCS417, pseudo1/3/5/6 |
| VillalobosEscobedo23 | https://doi.org/10.1371/journal.pgen.1010909 | (same) | HerbieS, Koxy, Putida, WCS417 |

### Table 1B. Linked publications not registered in feba.db

| pubId | URL (feba.db) | Published version (verified) | Linked strains (orgId) |
|---|---|---|---|
| **Georgoulis21** | (not registered) | mSystems 10.1128/msystems.01229-21 (PMC8562481; first author Georgoulis SJ, 2021) | **RalstoniaGMI1000, RalstoniaBSBF1503, RalstoniaPSI07** |

### Table 1C. Strains with no linked publication (not found by literature search either)
| orgId | Scientific name | Status |
|---|---|---|
| RalstoniaUW163 | R. solanacearum UW163 | **No linked publication** — Georgoulis21 includes only the other 3 strains; search for a UW163 RB-TnSeq paper found nothing |
| Methanococcus_JJ | Methanococcus maripaludis JJ | pubId empty + **archaeon → excluded** |
| Methanococcus_S2 | M. maripaludis S2 | 〃 |

---

## Table 2. Recovery (construction) medium and conditions per strain

- **Evidence:** recovery medium and temperature = original papers (the `Media` and selection-temperature columns of Price18 supplementary **Table S20_Mutagenesis** → 32 strains; Methods of individual papers → 14 strains). Oxygen = observed in feba.db `Experiment` (Price18 does not state it in the paper, so feba is the only source).
- **Evidence for vessel and shaking:** S20 has **no** vessel or shaking columns. Therefore, only what the Methods of each row's `pubId` paper state for the **recovery (thaw → mid-log phase) step** was transcribed.
  - `Price18`·`Price18 (Wetmore15 library)` 32 strains → common protocol in the main text: "an aliquot of the mutant library was thawed and inoculated into 25 ml rich medium with kanamycin or erythromycin and grown to mid-log phase **in a flask**" → **flask (25 mL)**; shaking rpm not stated.
  - `Georgoulis21` 3 strains → "three 1-ml aliquots of each transposon library ... revived in separate **100-ml flasks** of CPG medium with 12.5 mg/ml kanamycin **with shaking** incubation at 28°C for 16 to 20 h" → **100 mL flask + shaking yes** (rpm not stated. Routine subculturing in the same paper uses 200 rpm, but no rpm is given for the library recovery step).
  - `Helmann19`·`Helmann19b` → "inoculated into 25 mL fresh KB ... grown for ~7 h at 28°C **with shaking**" → vessel not stated, **shaking yes** (rpm not stated).
  - The rest (`Price18b`·`Liu21`·`Helmann21`·`Trotter21`·`McCausland21`·`Thompson18`·`Melnyk15`·UW163) → no vessel or shaking information at all for the recovery step, so **not stated**. Where only the culture volume is stated, the volume is kept in parentheses (e.g. not stated (50 mL)). The tubes, microplates and 200 rpm reported by Trotter21, Thompson18, etc. are conditions of the **fitness-assay step**, so they belong in Table 3, not Table 2.
- **Medium name normalization:** MOLS4=MoYLS4, marine broth=marine_broth_2216, ALP=ALP_rich.
> **2 archaeal strains excluded**: Methanococcus_JJ·Methanococcus_S2 (M. maripaludis) — archaea + no linked publication/essentiality source (Table 1C).

(46 bacterial strains, orgID in alphabetical order)

| orgId | Library | Recovery medium | Temp (°C) | Vessel | Shaking | Oxygen | pubId | Source text (summary) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| acidovorax_3H11 | `acidovorax_3H11_ML3a` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| ANA3 | `ANA3_ML1` | LB | 37 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=37 |
| azobra | `AzoBra_ML2a` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| BFirm | `BFirm_ML3` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| Btheta | `Btheta_ML6` | BHIS(+erythromycin) | 37 | not stated (50 mL) | not stated | **Anaerobic** | Liu21 | "thawed an aliquot ... into 50 mL BHIS supplemented with 10 μg/mL erythromycin ... incubating anaerobically at 37°C until mid-log (16 h) ... 'Time0'" |
| Burk376 | `Burk376_ML3` | LB(+kan 100) | 30 | not stated (50 mL) | not stated | Aerobic | Price18b | "at 30°C ... subinoculated into 50 ml LB with 100 µg/ml kanamycin ... recovered to mid-log (OD600~1.0)"; "cultured in LB during recovery from glycerol stocks" |
| Caulo | `Caulo_ML2` | PYE | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=PYE, selection temp=30 |
| Cola | `Cola_ML5` | marine broth (2216) | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=marine broth, selection temp=30 |
| Cup4G11 | `cupriavidus_4G11_ML11` | R2A | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=R2A, selection temp=30 |
| Dda3937 | `Dda3937` | LB(+kan) | 28 | not stated (25 mL) | not stated | Aerobic | Helmann21 | "recovered in 25 ml LB containing kanamycin at 28°C until OD600=0.5–0.7 ... two 1 ml pellets frozen as time0" |
| Ddia6719 | `Ddia6719` | LB(+kan) | 28 | not stated (25 mL) | not stated | Aerobic | Helmann21 | 〃 (same protocol for all 3 strains) |
| DdiaME23 | `DdiaME23` | LB(+kan) | 28 | not stated (25 mL) | not stated | Aerobic | Helmann21 | 〃 |
| Dino | `Dino_ML3` | marine broth (2216) | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=marine broth, selection temp=30 |
| DvH | `DvH_JW710` | MOYLS4 (lactate 60mM+sulfate 30mM+yeast extract 1g/L)(+G418 400µg/ml) | 30 | not stated (inside anaerobic chamber) | not stated | **Anaerobic** | Trotter21 | "inoculated in MOYLS4 medium supplemented with G418 (400 µg/mL) ... grow anaerobically at 30°C until mid-log ... 'Time0'" (anaerobic chamber 2%H₂/5%CO₂/93%N₂) |
| Dyella79 | `Dyella79_ML3` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| HerbieS | `HerbieS_ML4` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| Kang | `Kang_ML4` | marine broth (2216) | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=marine broth, selection temp=30 |
| Keio | `KEIO_ML9` | LB | 37 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=37 |
| Korea | `korea_ML2` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| Koxy | `Koxy_ML2` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| Magneto | `Magneto_ML2` | MG (+30 µM ferric malate, +kan 7µg/ml) | 30 | not stated (250 mL) | not stated | **Microaerobic (<10% O₂)** | McCausland21 | "standard laboratory conditions: MG medium supplemented with 30 µM iron and grown microaerobically ... one aliquot of ML2 thawed and inoculated into 250 ml MG with 7 µg/ml kanamycin ... time zero" |
| Marino | `Marino_ML2` | marine broth (2216) | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=marine broth, selection temp=30 |
| Miya | `Miya_ML3` | MoYLS4 (S20 notation MOLS4) | 30 | flask (25 mL) | not stated | **Anaerobic** | Price18 | S20: Media=MOLS4, selection temp=30 (note 7 = direct selection in 1 L MOLS4 + G418 liquid) |
| MR1 | `MR1_ML3` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 (Wetmore15 library) | S20: Media=LB, selection temp=30 (ref PMID25968644=Wetmore15) |
| Pedo557 | `Pedo557_ML3` | R2A | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=R2A, selection temp=30 |
| Phaeo | `Phaeo_ML1` | marine broth (2216) | 25 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=marine broth, selection temp=25 |
| Ponti | `Ponti_ML7` | marine broth (2216) | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=marine broth, selection temp=30 |
| PS | `PS_ML1` | ALP (ALP_rich) | 37 | not stated | not stated | Aerobic | Melnyk15 | S20: Media=ALP, selection temp=37 (electroporation, ref PMID25968643=Melnyk15) |
| pseudo1_N1B4 | `pseudo1_N1B4_ML1` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| pseudo3_N2E3 | `pseudo3_N2E3_ML2` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| pseudo5_N2C3_1 | `pseudo5_N2-C3_1_ML2` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| pseudo6_N2E2 | `pseudo6_N2E2_ML5` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| pseudo13_GW456_L13 | `pseudo13_ML2` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| psRCH2 | `psRCH2_ML7` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 (Wetmore15 library) | S20: Media=LB, selection temp=30 |
| Putida | `Putida_ML5_JBEI` | LB(+kan) | 30 | not stated (25 mL) | not stated | Aerobic | Thompson18 | "JBEI-1 thawed on ice and diluted in 25 ml LB supplemented with kanamycin ... time zero; monoculture at 30°C" (original construction paper = reference 16) |
| PV4 | `PV4_ML1` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| RalstoniaBSBF1503 | `RalstoniaBSBF1503_ML2` | CPG | 28 | 100 mL flask | yes (rpm not stated) | Aerobic | Georgoulis21 | "Ralstonia routinely grown in CPG rich medium (per L: 1g Casamino Acids, 10g Bacto peptone, 5g glucose, 1g yeast extract) at 28°C, shaking 200 rpm ... libraries revived in CPG" |
| RalstoniaGMI1000 | `RalstoniaGMI1000_ML2` | CPG (+kan 25µg/ml) | 28 | 100 mL flask | yes (rpm not stated) | Aerobic | Georgoulis21 | 〃 (same protocol for all 3 strains) |
| RalstoniaPSI07 | `RalstoniaPSI07_ML2` | CPG | 28 | 100 mL flask | yes (rpm not stated) | Aerobic | Georgoulis21 | 〃 |
| RalstoniaUW163 | `RalstoniaUW163_ML2` | CPG | 28 | not stated | not stated | Aerobic | (no dedicated RB-TnSeq paper; Lowe-Power lab CPG protocol) | "CPG broth ... 28°C" — standard medium for R. solanacearum (phylotype II, including UW163) (lowepowerlab.ucdavis.edu) |
| SB2B | `SB2B_ML5` | LB | 37 | flask (25 mL) | not stated | Aerobic | Price18 (Wetmore15 library) | S20: Media=LB, selection temp=37 (ref PMID25968644) |
| Smeli | `Smeli_ML6` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |
| SynE | `SynE_ML6` | BG11 (conjugation only in LB+DAP) | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=BG11, selection temp=30; note 5 = only SynE was conjugated in LB |
| SyringaeB728a | `SyringaeB728a_ML2` | KB (King's B)(+kan 100) | 28 | not stated (25 mL) | yes (rpm not stated) | Aerobic | Helmann19 | "inoculated from −80°C into 25 mL fresh KB with 100 µg/mL kanamycin, grown ~7 h at 28°C ... Time0; library generated in LB, inoculum·Time0 in KB" |
| SyringaeB728a_mexBdelta | `SyringaeB728a_mexBdelta_ML3` | KB (+kan 100) | 28 | not stated (25 mL) | yes (rpm not stated) | Aerobic | Helmann19b | "inoculated into 25 ml fresh KB with 100 µg/ml kanamycin and grown ~7 h at 28°C" |
| WCS417 | `WCS417_ML3` | LB | 30 | flask (25 mL) | not stated | Aerobic | Price18 | S20: Media=LB, selection temp=30 |

---

## Table 3. Standard unstressed condition per strain

This table compiles each strain's **standard unstressed (control) growth condition** from feba.db `Experiment`. Row order and orgIDs follow Table 2; cells where the Table 3 value differs from Table 2 (based on the literature/original sources) carry *(recovery: value)* alongside.

### Selection criteria

**0. Candidate set** — one **(library, medium, temperature, vessel, shaking) combination** forms one experiment candidate set. Only experiments identical in all five are grouped under one label.
- `mutantLibrary` is filled in for every experiment, but `media`·`vessel`·`temperature`·`shaking` are blank in 0.9%·0.7%·2.9%·3.7% of experiments, respectively (temperature and shaking are not recorded for the 3 Dickeya strains, and shaking is not recorded for Magneto). A blank is treated as a single value.
- Splitting is done **per experiment (expName)**; differing only in `pubId`, experimenter, or measurement date is not a reason to split.

**1. Exclusion of stress experiments** — the `expGroup` name, `condition_1~4`, and `expDesc` are examined **together**. Looking at only one of them is not sufficient.
- `condition_1~4` must be blank, and `expDesc` must contain no indication of a stress (pH value, drug concentration, starvation/survival, etc.).
- Experiments are excluded if `expGroup` belongs to the `stress`·`pH`·`survival`·`starvation`·`mouse`·`plant`·`in planta`·`xylem sap`·`motility`·`magnetic pulldown`·`metal limitation`·`vitamin`·`supernatant`·`reactor`·`anaerobic`·`fermentative growth` families. In some cases the stress is revealed only by the name even though the condition columns are blank (Btheta `set5IT061`·`set5IT062` are 48-hour aerobic exposures in the `survival` group).
- Conversely, an experiment may actually be unstressed even if its name is `temperature`·`nutrient`·`carbon source`, so experiments are not discarded by name alone (Smeli `set2IT069`·`set2IT070` have `expGroup`=`temperature` but are unstressed LB 30°C).
- **Exception** — selection markers originally contained in the recovery medium (kanamycin·G418·erythromycin) and the default carbon source of a defined medium (D-glucose) are not regarded as stresses.

**2. Eligibility** — a candidate must satisfy all three of the following to be considered for selection.

**(i) Exclusion of quality-outlier Time0 sets (†)** — within one candidate set, if any of `gMed`·`mad12`·`maxFit` deviates **3-fold or more relative to another Time0 set**, the worse set is excluded. "Worse" means lower `gMed` and higher `mad12`·`maxFit`, and a set that is worse in two or more of the three metrics is excluded. A low `gMed` means fewer BarSeq counts per gene, so the variance of the log ratios is larger; a high `mad12` means that insertions in the first and second halves of the same gene give different answers; and a high `maxFit` means there is a large gain that cannot arise in unstressed culture (a selection bottleneck or clonal sweep).

**(ii) $r_{cai} < 0.1$** — $r_{cai}$ = Spearman(codon adaptation index, `fit_mean`), computed over significant genes with $\lvert comb\_t \rvert \ge 2$. The codon adaptation index is a **sequence-only proxy for expression level**, calculated with reference to the codon usage of that strain's ribosomal proteins; because highly expressed genes are generally important, this correlation is normally negative. The threshold is set at 0.1 because using the sign alone rejects too many. The calculation is restricted to significant genes rather than all genes because the factor by which genes indistinguishable from 0 dilute the correlation differs by up to about 11-fold between strains.

| Verdict | Criterion | Action |
|---|---|---|
| Pass | $r_{cai} < 0$ | Adopted |
| Caution | $0 \le r_{cai} < 0.1$ | Adopted (flag retained) |
| **Inverted** | $r_{cai} \ge 0.1$ | Rejected |

**(iii) At least 100 significant genes** — the minimum among adopted experiment sets is 155 (Dino).

**3. Selection** — among the eligible candidates, the keys below are **compared in order, starting from the first**, and the single top-ranked candidate is chosen.

| Order | Key | Reason |
|---|---|---|
| ① | Medium and library **both** match Table 2 | The recovery condition reported in the literature is used as is |
| ② | Medium matches Table 2 or is a **rich medium** | The medium determines which genes are important, and 44 of the 46 strains are on rich media, so this keeps the labels homogeneous |
| ③ | Library matches Table 2 | The recovery step in Table 2 is precisely the step that produced that library's Time0 |
| ④ | Smaller difference from the recovery temperature, $\lvert \Delta T \rvert$ | Temperature is a real change of condition |
| ⑤ | Gentler shaking | Vigorous shaking was what flipped the phenotype axis. `orbital` records only the **mode** (orbital shaking), not an rpm, so it cannot be compared in magnitude with a numeric rpm; it is therefore placed as the gentlest |
| ⑥ | More significant genes | Statistical power |

② precedes ③ because changing the medium changes **the phenotype that `fit_mean` measures itself** (in a defined medium, amino-acid biosynthesis genes become important), whereas changing the library means **measuring the same phenotype with a different set of insertions**. Since each of the 389 `timeZeroSet`s belongs to a single library, fitness ratios are self-contained within whichever library is used, and a library mismatch with Table 2 is a mismatch of **documentation**, not of calculation. This does not mean, however, that substituting the library is harmless — the `fit_mean` agreement between candidate-set pairs that differ only in library is 0.16–0.45, comparable to substituting the medium. If ③ were removed, Btheta would miss its recovery library.

If ② were removed, 4 strains would move to a defined medium. The relative order of ① and ③ does not change the result.

### Time0-related columns

`Library` = `Experiment.mutantLibrary`. `No. Time0` = the number of distinct `Experiment.timeZeroSet` values referenced by those experiments. `Experiments per Time0` = `{timeZeroSet:[expName, ...], ...}`.

(46 bacterial strains, same order as Table 2)

| orgId | Library | Medium | Temp (°C) | Vessel | Shaking | Oxygen | Stress condition | pubId | No. experiments | No. Time0 | Experiments per Time0 | $r_{cai}$ (verdict) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| acidovorax_3H11 | `acidovorax_3H11_ML3a` | LB | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none | Price18 | 2 | 1 | {5/11/2015 acidovorax_3H11_ML3a_set3:[set3IT059, set3IT060]} | -0.157 pass |
| ANA3 ★ | `ANA3_ML1` | LB | 28 *(recovery: 37)* | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18 | 1 | 1 | {1/14/2015 ANA3_ML1_set4:[set4IT075]} | -0.111 pass |
| azobra ★ | `AzoBra_ML2a` | LB | 25 *(recovery: 30)* | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18 | 2 | 1 | {3/5/2015 AzoBra_ML2a_set1:[set1IT052, set1IT060]} | -0.295 pass |
| BFirm | `BFirm_ML3` | LB | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none | Price18 | 1 | 1 | {6/2/2015 BFirm_ML3_set2:[set2IT067]} | -0.175 pass |
| Btheta ‡§ | `Btheta_ML6` | Varel_Bryant_medium *(recovery: BHIS(+erythromycin))* | 37 | `tube` | `0 rpm` | **Anaerobic** | D-Glucose | Liu21 | 4 | 2 | {11-Oct-17 Btheta_ML6_set4:[set4IT047, set4IT050], 8-May-17 Btheta_ML6_set1:[set1IT071, set1IT072]} | -0.048 pass |
| Burk376 ‡§ | `Burk376_ML3` | RCH2_defined_noCarbon *(recovery: LB(+kan 100))* | 30 | `24-well transparent microplate; Multitron` | `700 rpm` | Aerobic | D-Glucose | Price18b | 4 | 2 | {1-May-17 Burk376_ML3_set1:[set1IT043, set1IT044], 5/2/18 Burk376_ML3_set4:[set4IT079, set4IT088]} | -0.188 pass |
| Caulo | `Caulo_ML2` | PYE | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none |  | 3 | 1 | {8/12/19 Caulo_ML2_set8:[set8IT011, set8IT023, set8IT035]} | +0.041 caution |
| Cola | `Cola_ML5` | marine_broth_2216 | 30 | `24-well transparent microplate; Multitron` *(recovery: flask)* | `700 rpm` | Aerobic | none | Price18 | 9 | 2 | {29-Mar-17 Cola_ML5_set1:[set1IT035, set1IT048, set1IT079], 8-May-17 Cola_ML5_set2:[set2IT011, set2IT024, set2IT034, set2IT045, set2IT058, set2IT068]} | -0.191 pass |
| Cup4G11 | `cupriavidus_4G11_ML11` | R2A | 30 | `24 deep-well microplate; Multitron` *(recovery: flask)* | `750 rpm` | Aerobic | none | Price18 | 3 | 2 | {5/13/2015 cupriavidus_4G11_ML11_set4:[set4IT069], 5/13/2015 cupriavidus_4G11_ML11_set5:[set5IT004, set5IT012]} | -0.246 pass |
| Dda3937 ★ | `Dda3937` | LB | not stated *(recovery: 28)* | `24 well plate` | not stated | **not stated** *(recovery: Aerobic)* | kanamycin | /Helmann21 | 9 | 3 | {1/26/21 Dda3937:[setATCH191, setATCH192], 6/5/21 Dda3937_set9:[set9TCH417, set9TCH418, set9TCH419, set9TCH420], 8/2/20 Dda3937:[setATCH110, setATCH111, setATCH112]} | -0.041 pass |
| Ddia6719 ★ | `Ddia6719` | LB | not stated *(recovery: 28)* | `24 well plate` | not stated | **not stated** *(recovery: Aerobic)* | kanamycin | Helmann21 | 4 | 1 | {1/26/21 Ddia6719:[setATCH203, setATCH204, setATCH206, setATCH207]} | -0.338 pass |
| DdiaME23 ★ | `DdiaME23` | LB | not stated *(recovery: 28)* | `24 well plate` | not stated | **not stated** *(recovery: Aerobic)* | kanamycin | /Helmann21 | 8 | 3 | {1/26/21 DdiaME23:[setATCH240, setATCH241, setATCH242], 6/5/21 DdiaME23_set10:[set10TCH436, set10TCH438], 8/2/20 DdiaME23:[setATCH113, setATCH114, setATCH115]} | -0.302 pass |
| Dino † | `Dino_ML3` | marine_broth_2216 | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none | Price18 | 2 | 1 | {1/26/2015 Dino_ML3_set7:[set7IT086, set7IT087]} | -0.054 pass |
| DvH | `DvH_JW710` | MoYLS4 | 30 | `24 deep-well microplate` | `0 rpm` | **Anaerobic** | none | Trotter21 | 12 | 4 | {11/10/2016 DvH_JW710_set4:[set4IT040, set4IT041, set4IT042], 4/18/2016 DvH_JW710_set2:[set2S268, set2S269, set2S270], 4/6/2016 DvH_JW710_set1:[set1S209, set1S210, set1S211], 5/5/2016 DvH_JW710_set3:[set3S344, set3S345, set3S346]} | -0.009 pass |
| Dyella79 | `Dyella79_ML3` | LB | 30 | `24-well transparent microplate; Multitron` *(recovery: flask)* | `700 rpm` | Aerobic | none | Price18 | 4 | 1 | {17-May-17 Dyella79_ML3_set3:[set3IT040, set3IT051, set3IT060, set3IT067]} | -0.300 pass |
| HerbieS ‡‖ | `HerbieS_ML4a` *(recovery: `HerbieS_ML4`)* | R2A *(recovery: LB)* | 30 | `24-well transparent microplate; Multitron` *(recovery: flask)* | `700 rpm` | Aerobic | none | VillalobosEscobedo23 | 6 | 2 | {1-Feb-22 HerbieS_ML4_set13:[set13IT052, set13IT053, set13IT054], 17-Aug-21 HerbieS_ML4_set12:[set12IT038, set12IT039, set12IT040]} | -0.078 pass |
| Kang † | `Kang_ML4` | marine_broth_2216 | 30 | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18 | 4 | 1 | {3/12/2014 Kang_ML4_set2:[set2H11, set2H20, set2H43, set2H48]} | -0.067 pass |
| Keio ★ | `Keio_ML9` | LB | 28 *(recovery: 37)* | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18/Wetmore15 | 2 | 1 | {6/30/2014 Keio_ML9_set2:[set2IT045, set2IT069]} | -0.125 pass |
| Korea | `korea_ML2` | LB | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none | Price18 | 2 | 1 | {1/26/2015 korea_ML2_set5:[set5IT075, set5IT076]} | -0.052 pass |
| Koxy | `Koxy_ML2` | LB | 30 | `24-well transparent microplate; Multitron` *(recovery: flask)* | `700 rpm` | Aerobic | none | Price18 | 5 | 2 | {27-Apr-17 Koxy_ML2_set6:[set6IT088], 27-Apr-17 Koxy_ML2_set7:[set7IT016, set7IT054, set7IT071, set7IT078]} | -0.111 pass |
| Magneto | `Magneto_ML2` | MG minimal media | 30 | `tube` | not stated | **Microaerobic** *(recovery: **microaerobic (<10% O₂)**)* | none | McCausland21 | 1 | 1 | {27-Sep-18 Magneto_Ml2_set1:[set1IT053]} | -0.185 pass |
| Marino | `Marino_ML2` | marine_broth_2216 | 30 | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18 | 11 | 4 | {10/23/2013 Marino_ML2_set4:[set4H10, set4H18, set4H30], 9/10/2013 Marino_ML2_set3:[set3H25, set3H7], 9/11/2013 Marino_ML2_set3:[set3H30, set3H31, set3H39], 9/9/2013 Marino_ML2_set2:[set2H22, set2H31, set2H41]} | -0.151 pass |
| Miya | `Miya_ML3` | MoYLS4 | 30 | `24 deep-well microplate` *(recovery: flask)* | `0 rpm` | **Anaerobic** | none | Price18 | 3 | 1 | {11/23/2014 Miya_ML3_set2:[set2IT012, set2IT016, set2IT025]} | -0.163 pass |
| MR1 | `MR1_ML3` | LB | 30 | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Wetmore15 | 12 | 6 | {10/22/2013 MR1_ML3_set6:[set6H11, set6H18, set6H25, set6H7], 8/18/2013 MR1_ML3_set1:[set1H15], 8/18/2013 MR1_ML3_set2:[set2H29, set2H44], 8/19/2013 MR1_ML3_set4:[set4H43, set4H48], 8/19/2013 MR1_ML3_set5:[set5H12], 8/21/2013 MR1_ML3_set5:[set5H44, set5H47]} | -0.166 pass |
| Pedo557 | `Pedo557_ML3` | R2A | 30 | `24-well transparent microplate; Multitron` *(recovery: flask)* | `700 rpm` | Aerobic | none | Price18 | 9 | 2 | {17-May-17 Pedo557_ML3_set1:[set1IT019, set1IT034, set1IT053, set1IT068, set1IT095], 17-May-17 Pedo557_ML3_set2:[set2IT008, set2IT022, set2IT032, set2IT041]} | -0.225 pass |
| Phaeo | `Phaeo_ML1` | marine_broth_2216 | 25 | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18/Wetmore15 | 15 | 7 | {1/20/2015 Phaeo_ML1_set8:[set8IT088, set8IT096], 8/1/2013 Phaeo_ML1_set2:[set2H31, set2H5], 8/5/2013 Phaeo_ML1_set2:[set2H43], 8/5/2013 Phaeo_ML1_set3:[set3H7], 8/6/2013 Phaeo_ML1_set4:[set4H18], marchapr14:[set7IT053, set7IT062, set7IT073, set7IT085], marchapr14 Phaeo_ML1_set6:[set6IT009, set6IT016, set6IT030, set6IT041]} | -0.099 pass |
| Ponti | `Ponti_ML7` | marine_broth_2216 | 30 | `24-well transparent microplate; Multitron` *(recovery: flask)* | `700 rpm` | Aerobic | none | Price18 | 10 | 3 | {29-Mar-17 Ponti_ML7_set1:[set1IT035, set1IT045, set1IT063, set1IT073], 8-May-17 Ponti_ML7_set2:[set2IT084, set2IT096], 8-May-17 Ponti_ML7_set3:[set3IT015, set3IT028, set3IT040, set3IT055]} | -0.274 pass |
| PS | `PS_ML1` | ALP_rich | 37 | `tube` | `200 rpm` | Aerobic | none | Melnyk15 | 3 | 1 | {8/1/2013 PS_ML1_set1:[set1H13, set1H15, set1H16]} | -0.226 pass |
| pseudo1_N1B4 | `pseudo1_N1B4_ML1` | LB | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none | Price18 | 1 | 1 | {10/18/2014_or_07/22/16 pseudo1_N1B4_ML1_set5:[set5IT018]} | -0.162 pass |
| pseudo3_N2E3 | `pseudo3_N2E3_ML2` | LB | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none | Price18 | 2 | 1 | {10/28/2015 pseudo3_N2E3_ML2_set6:[set6IT053, set6IT054]} | -0.325 pass |
| pseudo5_N2C3_1 | `pseudo5_N2-C3_1_ML2` | LB | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none | Price18 | 2 | 1 | {10/28/2014 pseudo5_N2-C3_1_ML2_set6:[set6IT074, set6IT075]} | -0.174 pass |
| pseudo6_N2E2 | `pseudo6_N2E2_ML5` | LB | 30 | `24 deep-well microplate; Multitron` *(recovery: flask)* | `750 rpm` | Aerobic | none | Price18 | 4 | 1 | {3/14/2015 pseudo6_N2E2_ML5_set3:[set3IT072, set3IT076, set3IT082, set3IT085]} | -0.115 pass |
| pseudo13_GW456_L13 ★ | `pseudo13_ML2` | LB | 25 *(recovery: 30)* | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18 | 1 | 1 | {3/11/2015 pseudo13_ML2_set2:[set2IT087]} | -0.299 pass |
| psRCH2 | `psRCH2_ML7` | LB | 30 | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18/Wetmore15 | 6 | 4 | {6/24/2013 psRCH2_ML7_set3:[set3H43], 6/24/2013 psRCH2_ML7_set4:[set4H19, set4H6], 6/26/2013 psRCH2_ML7_set4:[set4H37], 6/28/2013 psRCH2_ML7_set5:[set5H14, set5H5]} | -0.143 pass |
| Putida ‡‖ | `Putida_ML5` *(recovery: `Putida_ML5_JBEI`)* | R2A *(recovery: LB(+kan))* | 30 | `24-well transparent microplate; Multitron` | `700 rpm` | Aerobic | none | VillalobosEscobedo23 | 6 | 2 | {1-Feb-22 Putida_ML5_set36:[set36IT074, set36IT075, set36IT076], 28-Jun-22 Putida_ML5_set38:[set38IT007, set38IT010, set38IT011]} | -0.194 pass |
| PV4 | `PV4_ML1` | LB | 30 | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Price18 | 4 | 3 | {12/11/2013 PV4_ML1_set1:[set1H25], 12/11/2013 PV4_ML1_set3:[set3IT055], 2/20/2014 PV4_ML1_set1:[set1H1, set1H6]} | +0.063 caution |
| RalstoniaBSBF1503 | `RalstoniaBSBF1503_ML2` | CPG | 28 | `tube` *(recovery: 100 mL flask)* | `200 rpm` | Aerobic | none |  | 3 | 2 | {8-Feb-19 RalstoniaBSBF1503_ML2_set1:[set1IT076, set1IT080], 8-Feb-19 RalstoniaBSBF1503_ML2_set3:[set3IT039]} | -0.252 pass |
| RalstoniaGMI1000 | `RalstoniaGMI1000_ML2` | CPG | 28 | `tube` *(recovery: 100 mL flask)* | `200 rpm` | Aerobic | none |  | 2 | 2 | {8-Feb-19 RalstoniaGMI1000_ML2_set2:[set2IT062], 8-Feb-19 RalstoniaGMI1000_ML2_set4:[set4IT033]} | -0.119 pass |
| RalstoniaPSI07 | `RalstoniaPSI07_ML2` | CPG | 28 | `tube` *(recovery: 100 mL flask)* | `200 rpm` | Aerobic | none |  | 2 | 1 | {8-Feb-19 RalstoniaPSI07_ML2_set1:[set1IT085, set1IT089]} | -0.139 pass |
| RalstoniaUW163 | `RalstoniaUW163_ML2` | CPG | 28 | `tube` | `200 rpm` | Aerobic | none |  | 3 | 2 | {8-Feb-19 RalstoniaUW163_ML2_set1:[set1IT067, set1IT071], 8-Feb-19 RalstoniaUW163_ML2_set2:[set2IT036]} | -0.233 pass |
| SB2B | `SB2B_ML5` | LB | 37 | `48 well microplate; Tecan Infinite F200` *(recovery: flask)* | `orbital` | Aerobic | none | Wetmore15 | 14 | 5 | {10/15/2013 MR1_SB2B_combined_set:[setAH41], 3/25/2013 SB2B_ML5_set1:[set1H42, set1H43, set1H44], 3/25/2013 SB2B_ML5_set2:[set2H12, set2H6, set2H7], 7/26/2013 SB2B_ML5_set4:[set4H10, set4H20, set4H34, set4H42], 7/30/2013 SB2B_ML5_set5:[set5H11, set5H23, set5H7]} | -0.194 pass |
| Smeli | `Smeli_ML6` | LB | 30 | `tube` *(recovery: flask)* | `200 rpm` | Aerobic | none | Price18 | 2 | 1 | {6/2/2015 Smeli_ML6_set2:[set2IT069, set2IT070]} | -0.020 pass |
| SynE | `SynE_ML6` | BG11 | 30 | `12 well microplate` *(recovery: flask)* | `0 rpm` | Aerobic | none | Price18 | 18 | 4 | {10/2/2015 SynE_ML6_set7:[set7IT032, set7IT042, set7IT047], 10/2/2015 SynE_ML6_set8:[set8IT052, set8IT057, set8IT062, set8IT066, set8IT069, set8IT080, set8IT085, set8IT090], 6/20/2015 SynE_ML6_set3:[set3IT004, set3IT016, set3IT031], 6/20/2015 SynE_ML6_set4:[set4IT050, set4IT062, set4IT077, set4IT089]} | -0.106 pass |
| SyringaeB728a | `SyringaeB728a_ML2` | KB | 28 | `24 well microplate` | `250 rpm` | Aerobic | none | Helmann19 | 2 | 2 | {12/4/17 SyringaeB728a_ML2_set1:[set1IT005], 6/8/18 SyringaeB728a_ML2_set4:[set4IT013]} | -0.048 pass |
| SyringaeB728a_mexBdelta | `SyringaeB728a_mexBdelta_ML3` | KB | 28 | `24 well microplate` | `250 rpm` | Aerobic | none | Helmann19b | 2 | 2 | {12/4/17 SyringaeB728a_mexBdelta_ML3_set1:[set1IT006], 6/8/18 SyringaeB728a_mexBdelta_ML3_set4:[set4IT014]} | -0.028 pass |
| WCS417 ‡ | `fluoroDangl_ML3` | R2A *(recovery: LB)* | 30 | `24-well transparent microplate; Multitron` *(recovery: flask)* | `700 rpm` | Aerobic | none | VillalobosEscobedo23 | 1 | 1 | {17-Aug-21 fluoroDangl_ML3_set18:[set18IT044]} | +0.015 caution |

**Symbols** — mark strains that did not satisfy selection keys ①–④ or that had experiments excluded under eligibility 2-(i).
- **‡** = strains for which no eligible candidate existed on the recovery medium, so a **different medium** was chosen — Btheta·Burk376·HerbieS·Putida·WCS417.
- **§** = strains for which no candidate existed on a rich medium, so a **defined medium** was chosen — Btheta·Burk376.
- **‖** = strains for which no candidate existed with the recovery library, so a **different library** was chosen — HerbieS·Putida.
- **★** = 7 strains for which no candidate existed at the recovery temperature, or the temperature was not recorded, so a **different growth temperature** was chosen — ANA3·azobra·Dda3937·Ddia6719·DdiaME23·Keio·pseudo13_GW456_L13.
- **†** = strains with some experiments excluded as quality outliers under eligibility 2-(i) — Dino (1 experiment)·Kang (3 experiments).
- **(recovery: value)** = the cell differs from Table 2, so the Table 2 value is shown alongside. *Temperature*: Table 2 gives the construction/recovery temperature, Table 3 the growth temperature of the fitness assay. *Vessel*: shown only for strains whose recovery vessel is stated in the literature (Price18 series = flask, Georgoulis21 3 strains = 100 mL flask). *Shaking*: not shown, because no paper states the shaking speed of the recovery step.

**Caveats**
- **Library notation**: Table 2 takes library names from Price18 supplementary **Table S20_Mutagenesis** `Mutant library name` (32 strains), the paper main text (Putida=JBEI-1, DvH=JW710 background, Magneto=ML2), and feba.db `Experiment.mutantLibrary` (remaining 11 strains), in that order. `Keio` (S20 notation `KEIO_ML9`, differing only in letter case) and `WCS417` (renamed `fluoroDangl_ML3` in feba) are not shown alongside.
- **Vessel and shaking notation**: feba values are transcribed as is. `; Multitron`·`; Tecan Infinite F200` are the incubator/plate-reader models recorded in the same field. `tube`=test tube, `24 deep-well microplate`=24-well deep-well microplate, `orbital`=orbital shaking (rpm not recorded), `0 rpm`=static.
- **Medium name notation**: `marine_broth_2216`=Table 2 "marine broth (2216)", `ALP_rich`=Table 2 "ALP", `MoYLS4`=Table 2 "MOYLS4/MOLS4" — these are notation differences for the same medium, so they are not shown alongside.
- **Selection markers and supplements not shown**: the selection markers and supplements of the literature recovery media (Btheta erythromycin, DvH G418, SyringaeB728a kan 100, RalstoniaGMI1000 kan 25µg/ml, Magneto ferric malate + kan) are not recorded in feba `media`, so they do not appear in the Table 3 medium names. For Dickeya, kanamycin is recorded in the condition columns and is shown as is.
- **Single replicate**: ANA3·BFirm·Magneto·pseudo1_N1B4·pseudo13_GW456_L13·WCS417 have only one experiment each, so reproducibility across replicates cannot be checked.
- **$r_{cai}$ caution verdict**: Caulo·PV4·WCS417 have $r_{cai}$ of at least 0 and below 0.1, so they were adopted, but the flag is retained.
- **Caulo `pubId` blank**: the 3 adopted `Caulo_ML2_set8` experiments have an empty `pubId` and a different experimenter and measurement date (Justin Zik, 2019-08-12). The library is the same `Caulo_ML2`.


## Table 4. Essential gene information per strain

**Notation rule**: every numeric column uses the format **`all gene types (protein-coding)`** (the value in parentheses is type=1 protein-coding). Example: `4610 (4146)` = 4146 protein-coding genes out of 4610 genes of all types. The values inside and outside the parentheses are **each independently consistent**, so the arithmetic can be checked on either basis.

$\text{unknown} = \text{total} - \text{essential} - \text{measured}$ (on the all-types basis and the coding basis, separately)

- **`total`** = number of feba.db `Gene` entries.
- **`nonessential(measured)`** = number of unique genes with a fitness value in feba.db `GeneFitness`. Essential genes, by definition, have no fitness value. **Non-coding genes also receive a fitness value when insertions occur in them** (pseudogenes type=7, tRNA type=5, rRNA type=2, etc.), numbering **1–655 per strain (median 12)**. The strains where they make up a large share are **azobra 655 (13.6% of measured)** · Keio 204 (5.4%) · DdiaME23 150 (4.0%) · the 4 Ralstonia strains 127–131 (3.2–3.4%) · MR1 120 (3.2%), while Btheta and HerbieS have only 1. In such strains, reading the value outside the parentheses as "the number of nonessential protein-coding genes" overestimates it, so **use the value inside the parentheses (coding basis)**.
- **`essential`** = the sources below (author reports or Price18 `essential_proteins.tab`) — joining the essential lists of all 42 strains to feba.db `Gene.type` showed that **all were type=1 (protein-coding), with 0 non-coding or unmatched genes** (`essential_proteins.tab` **13,869/13,869** (all 32 strains) + all 10 files in `refer/`). In other words, "essential" is, by definition, assigned only to protein-coding genes. *However, for Cup4G11, the species name in `essential_proteins.tab` (`Cupriavidus basilensis 4G11`) differs from that in feba.db `Organism` (`FW507-4G11`), so it was confirmed by **manual mapping***
- **`unknown`** = **number of genes called neither essential nor nonessential** (fitness could not be obtained because of insufficient reads, or the gene is not in the essential list).

### Table 4A. 42 strains with essentiality data (orgId ascending)

The 32 strains whose `Source` is `essential_proteins.tab` belong to the Price18 series (bigfit, http://genomics.lbl.gov/supplemental/bigfit/); the 10 strains with `refer/...` use author-reported essential lists from 2018 onward.

| orgId | essential | nonessential(measured) | total | unknown | Source |
|---|---|---|---|---|---|
| acidovorax_3H11 | 504 (504) | 3935 (3924) | 5019 (4964) | 580 (536) | essential_proteins.tab |
| ANA3 | 307 (307) | 3668 (3643) | 4552 (4360) | 577 (410) | essential_proteins.tab |
| azobra | 340 (340) | 4833 (4178) | 6447 (5488) | 1274 (970) | essential_proteins.tab |
| BFirm | 409 (409) | 5428 (5422) | 7263 (7182) | 1426 (1351) | essential_proteins.tab |
| Btheta | 378 (378) | 4055 (4054) | 4902 (4816) | 469 (384) | refer/Liu21_Btheta_essential.xlsx (Cell Rep Table S1) |
| Caulo | 440 (440) | 3312 (3305) | 3943 (3886) | 191 (141) | essential_proteins.tab |
| Cola | 492 (492) | 3954 (3950) | 4684 (4625) | 238 (183) | essential_proteins.tab |
| Cup4G11 | 474 (474) | 6384 (6373) | 7449 (7358) | 591 (511) | essential_proteins.tab |
| Dda3937 | 374 (374) | 3705 (3614) | 4404 (4213) | 325 (225) | refer/Helmann21_Dda3937_essential.csv (GitHub tylerhelmann/dickeya-barseq-2021) |
| Ddia6719 | 426 (426) | 3528 (3433) | 4338 (4110) | 384 (251) | refer/Helmann21_Ddia6719_essential.csv |
| DdiaME23 | 426 (426) | 3761 (3611) | 4437 (4182) | 250 (145) | refer/Helmann21_DdiaME23_essential.csv |
| Dino | 535 (535) | 3187 (3185) | 4562 (4192) | 840 (472) | essential_proteins.tab |
| DvH | 399 (399, wt) | 2741 (2725) | 3668 (3403) | 528 (279) | refer/Trotter21_DvH_essential.xlsx (Front. Table S5, ess.wt) |
| Dyella79 | 371 (371) | 3625 (3611) | 4395 (4317) | 399 (335) | essential_proteins.tab |
| HerbieS | 472 (472) | 3879 (3878) | 4779 (4715) | 428 (365) | essential_proteins.tab |
| Kang | 399 (399) | 2003 (1989) | 2513 (2463) | 111 (75) | essential_proteins.tab |
| Keio | 324 (324) | 3789 (3585) | 4610 (4146) | 497 (237) | essential_proteins.tab |
| Korea | 443 (443) | 3393 (3360) | 4245 (4149) | 409 (346) | essential_proteins.tab |
| Koxy | 378 (378) | 4608 (4485) | 5586 (5309) | 600 (446) | essential_proteins.tab |
| Marino | 555 (555) | 3650 (3648) | 4470 (4410) | 265 (207) | essential_proteins.tab |
| Miya | 312 (312) | 2531 (2511) | 3601 (3180) | 758 (357) | essential_proteins.tab |
| MR1 | 397 (397) | 3782 (3662) | 4745 (4467) | 566 (408) | essential_proteins.tab |
| Pedo557 | 478 (478) | 4423 (4359) | 5089 (4964) | 188 (127) | essential_proteins.tab |
| Phaeo | 538 (538) | 3099 (3094) | 3944 (3875) | 307 (243) | essential_proteins.tab |
| Ponti | 472 (472) | 3685 (3643) | 4329 (4220) | 172 (105) | essential_proteins.tab |
| PS | 485 (485) | 2561 (2550) | 3504 (3436) | 458 (401) | essential_proteins.tab |
| pseudo1_N1B4 | 426 (426) | 4336 (4333) | 6058 (5972) | 1296 (1213) | essential_proteins.tab |
| pseudo3_N2E3 | 405 (405) | 5028 (5024) | 5854 (5766) | 421 (337) | essential_proteins.tab |
| pseudo5_N2C3_1 | 387 (387) | 5193 (5187) | 6082 (6000) | 502 (426) | essential_proteins.tab |
| pseudo6_N2E2 | 442 (442) | 5133 (5126) | 6175 (6094) | 600 (526) | essential_proteins.tab |
| pseudo13_GW456_L13 | 325 (325) | 4350 (4346) | 5243 (5152) | 568 (481) | essential_proteins.tab |
| psRCH2 | 430 (430) | 3349 (3345) | 4336 (4265) | 557 (490) | essential_proteins.tab |
| PV4 | 289 (289) | 3009 (2998) | 3993 (3859) | 695 (572) | essential_proteins.tab |
| RalstoniaBSBF1503 | 454 (454) | 3779 (3652) | 4897 (4659) | 664 (553) | refer/Georgoulis21_ST1_essential.xlsx (IBSBF1503 sheet) |
| RalstoniaGMI1000 | 456 (456) | 3980 (3853) | 5204 (4956) | 768 (647) | refer/Georgoulis21_ST1_essential.xlsx (GMI1000 sheet, Essentiality.pl=Price18 method) |
| RalstoniaPSI07 | 471 (471) | 4029 (3898) | 4968 (4735) | 468 (366) | refer/Georgoulis21_ST1_essential.xlsx (PSI07 sheet) |
| SB2B | 394 (394) | 3121 (3099) | 3785 (3645) | 270 (152) | essential_proteins.tab |
| Smeli | 559 (559) | 5133 (5130) | 6281 (6217) | 589 (528) | essential_proteins.tab |
| SynE | 614 (614) | 1899 (1898) | 2722 (2669) | 209 (157) | essential_proteins.tab |
| SyringaeB728a | 392 (392) | 4314 (4271) | 5216 (5089) | 510 (426) | refer/Helmann19_SyringaeB728a_essential.xlsx (PNAS Dataset S1, `ess` column) |
| SyringaeB728a_mexBdelta | 392 (392) | 4282 (4237) | 5216 (5089) | 542 (460) | shares the file above (same B728a genome) |
| WCS417 | 473 (473) | 4419 (4414) | 5580 (5506) | 688 (619) | essential_proteins.tab |

### Table 4B. 4 strains without essentiality data — total and measured only (observed in feba.db)

These strains are not included in the Price18 series (essential_proteins.tab), and no author-reported essential list is available in refer/, so **essential cannot be computed** (the raw Price18-style dens·normreads data are unavailable). Nonessential (measured) and total were computed from feba.db.

| orgId | essential | nonessential(measured) | total | Remarks |
|---|---|---|---|---|
| Burk376 | unknown | 5153 (5147) | 6609 (6555) | Paraburkholderia bryophila 376MFSha3.1; not included in bigfit essential_proteins.tab. **No essentiality data in the supplementary material** |
| Magneto | unknown (reported: 445) | 3221 (3155) | 4691 (4489) | Magnetospirillum magneticum AMB-1; **the McCausland21 main text reports 445 essential genes, but the supplementary material contains no gene list** |
| Putida | unknown | 4778 (4769) | 5661 (5563) | Pseudomonas putida KT2440; JBEI library (original construction paper not included in refer/). The nonessential count of 4778 agrees with the value reported in the literature |
| RalstoniaUW163 | unknown | 3795 (3667) | 5006 (4758) | R. solanacearum UW163; no linked publication (Table 1C) → no essentiality source |
