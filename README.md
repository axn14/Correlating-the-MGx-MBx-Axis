# CRC Gut Microbiome–Metabolome Integration Study

<img width="1193" height="667" alt="image" src="https://github.com/user-attachments/assets/a4d14ef4-da9f-4f43-a6b7-0ed70b243200" />
<img width="994" height="695" alt="image" src="https://github.com/user-attachments/assets/959f8f36-a9ac-46ec-8e6b-0cb6bef699b6" />


**Generalised multi-omics pipeline for identifying microbial and metabolic signatures across colorectal cancer progression.**

---

## Quick Start

### Option A — Browse pre-computed results (no installation required beyond Streamlit)

```bash
git clone https://github.com/axn14/Correlating-the-MGx-MBx-Axis.git
cd Correlating-the-MGx-MBx-Axis
pip install -r streamlit_app/requirements_st.txt
streamlit run streamlit_app/app.py
# Open http://localhost:8501
```

The results viewer loads pre-computed figures and tables from the YACHIDA-CRC-2019 cohort. No heavy dependencies required.

### Option B — Run the pipeline on your own data (Docker)

```bash
git clone https://github.com/axn14/Correlating-the-MGx-MBx-Axis.git
cd Correlating-the-MGx-MBx-Axis
docker-compose up --build
# Open http://localhost:8501 → click "Upload & Run" in the sidebar
```

Upload your own species abundance table, metabolomics table, and sample metadata. The pipeline runs steps 1–5 (Preprocessing → Associations → ML → MOFA+ → Network) fully inside the container. See [User Data Pipeline](#user-data-pipeline) for input format details and limitations.

> **⚠ Disclaimer:** This tool is provided for academic and research use only. It has been developed and validated exclusively on the YACHIDA-CRC-2019 cohort. Results on other datasets have not been clinically validated and should not be used for diagnostic or clinical decision-making purposes. See [Limitations](#limitations) for details on which analysis streams are unavailable for user-submitted data.

---

## Overview

This repository contains the complete analysis pipeline for an integrative microbiome–metabolome study of colorectal cancer (CRC). The pipeline combines shotgun metagenomic species profiles, untargeted LC-MS metabolomics, machine learning, constraint-based metabolic modelling, and genomic metabolic capacity inference (GutSMASH) to identify microbial sources of disease-associated metabolites across CRC progression stages.

### Primary Cohort
| Dataset | Samples | Stages | Data Types |
|---|---|---|---|
| **Yachida et al. 2019** (*Nature Medicine*) | 347 | Healthy, Stage 0, HS, Stage I/II, Stage III/IV, MP | Shotgun metagenomics + LC-MS metabolomics |

### Validation Cohort
| Dataset | Samples | Stages | Data Types |
|---|---|---|---|
| **Sinha et al. 2016** (*Science*) | — | Healthy vs CRC | Shotgun metagenomics |

### Planned LODO Cohorts
- **Kong et al. 2023** (*Nature Communications*) — early-onset CRC, shotgun + LC-MS
- **Sinha-Vogtmann shotgun** — shotgun replacement for Sinha 2016
- **Borenstein Lab CRC** — curated paired shotgun-metabolomic dataset

---

## Pipeline Architecture

```
Raw reads / mzML files
        │
        ▼
NB00 ── Upstream processing (fastp → MetaPhlAn4 → HUMAnN3 → mzML feature extraction)
NB00b ─ pyOpenMS LC-MS pipeline (FeatureFinderMetabo → adduct grouping → RT alignment)
        │
        ▼
NB01 ── Preprocessing & normalisation (CLR species, log₁₀ metabolites, QC, PCA)
        │
        ▼
NB02 ── Association analysis (DA, Spearman/partial correlations, stage-stratified, networks)
        │
        ▼
NB03 ── ML benchmarking (XGBoost, LightGBM, ElasticNet, SVR, RF · SHAP · LODO CV)
        │
        ▼
NB04 ── Multivariate trajectory (MOFA+, FDA, ordinal regression, MP recovery arc)
        │
        ▼
NB05 ── Mediation & pathway network (bipartite network, PERMANOVA, bootstrap ACME)
        │
        ▼
NB06 ── GutSMASH benchmarking (BGC genomic capacity, literature evidence, novel candidates)
        │
        ▼
NB07 ── Advanced evidence integration (multi-source scoring, cross-notebook synthesis)
        │
        ▼
NB08 ── Metabolite source attribution (generalised species-to-metabolite attribution)
        │
        ▼
NB09 ── Mechanistic flux (MICOM community FBA, AGORA103, excretion flux atlas)
```

**Holy Trinity integration:** SHAP (NB03, statistical) × GutSMASH BGC (NB06, genomic) × MICOM flux (NB09, mechanistic) combined in a single publication figure per metabolite.

---

## Repository Structure

```
├── Code/                                   # Same notebooks (copy for direct code browsing)
│   ├── 00_upstream_pipeline.ipynb          # SRA → fastp → MetaPhlAn4 → HUMAnN3 → mzML
│   ├── 00b_lcms_pyopenms.ipynb             # pyOpenMS LC-MS feature detection pipeline
│   ├── 01_preprocessing.ipynb              # CLR/log10 transform, QC, PCA, duplicates
│   ├── 02_association_analysis.ipynb       # DA, correlations, stage-stratified, networks
│   ├── 03_ml_benchmarking.ipynb            # 5-model CV, SHAP, LODO cross-validation
│   ├── 04_mofa_fda_trajectory.ipynb        # MOFA+, FDA, ordinal regression, MP trajectory
│   ├── 05_mediation_network.ipynb          # Bipartite network, PERMANOVA, bootstrap ACME
│   ├── 06_gutsmash_benchmarking.ipynb      # GutSMASH runner, BGC matching, literature tiers
│   ├── 07_advanced_evidence_integration.ipynb  # Multi-source evidence integration
│   ├── 08_general_metabolite_source_attribution.ipynb  # Generalised source attribution
│   └── 09_mechanistic_flux_micom.ipynb     # MICOM FBA, AGORA103 community models
├── pipeline/                               # Modular pipeline for user-submitted data
│   ├── data_intake.py                      # Format detection: MetaPhlAn, OTU, KEGG/HMDB
│   ├── step1_preprocess.py                 # CLR + log10 transforms, prevalence filters
│   ├── step2_associations.py               # DA, Spearman, partial correlations
│   ├── step3_ml.py                         # 5-model nested CV + SHAP
│   ├── step4_mofa.py                       # MOFA+ wrapper (mofapy2)
│   ├── step5_network.py                    # Partial-corr network + bootstrap mediation
│   ├── runner.py                           # Orchestrator with background-thread support
│   └── reference/                          # HMDB→KEGG and compound name lookup tables
├── streamlit_app/                          # Interactive results viewer
│   ├── app.py                              # Home page
│   ├── utils_st.py                         # Shared helpers
│   ├── requirements_st.txt                 # Minimal deps for viewer only
│   └── pages/
│       ├── 0_Upload_and_Run.py             # Upload UI + pipeline runner
│       ├── 1_Preprocessing.py
│       ├── 2_Associations.py
│       ├── 3_Machine_Learning.py
│       ├── 4_MOFA.py
│       ├── 5_Mediation_Network.py
│       ├── 6_Evidence_Integration.py
│       └── 7_MICOM_Flux.py
├── Thesis-Figures/pdf converted/           # Pre-computed publication figures (31 images)
├── Dockerfile                              # Full pipeline container (Python 3.12 + R)
├── docker-compose.yml                      # Volume-mounted run configuration
├── utils.py                                # Shared constants, transforms, DA, correlation
├── requirements.txt                        # pip-installable dependencies
├── environment.yml                         # conda environment (Python 3.12)
├── Data/                                   # Raw TSV files (not tracked; see Data section)
└── Results/                                # Output figures, tables, intermediates (not tracked)
```

---

## Methods Summary

### 1. Upstream Processing (NB00, NB00b) [BETA]
- **Quality control:** fastp (adapter trimming, quality filtering, paired-end)
- **Taxonomic profiling:** MetaPhlAn4 (mpa_vJan21_CHOCOPhlAnSGB_202103 database)
- **Functional profiling:** HUMAnN3 (pathway abundance + gene families)
- **LC-MS feature detection:** pyOpenMS `FeatureFinderMetabo` → `MetaboliteFeatureDeconvolution` (adduct grouping) → `MapAlignerPoseClustering` (RT alignment)
- **KEGG annotation:** REST API bulk lookup of m/z-matched compound IDs

### 2. Preprocessing (NB01)
- **Species:** centre-log ratio (CLR) transform with Aitchison pseudocount
- **Metabolites:** log₁₀ transform; KEGG ID–name resolution
- **QC filters:** prevalence ≥ 10%, variance bottom 10% removed per dataset
- **Metadata harmonisation:** 6-level CRC stage ordering (Healthy → MP); 3-group map (Healthy / Early_CRC / Advanced_CRC)
- **Duplicate detection:** ERAWIJANTARI–YACHIDA shared sample removal

### 3. Association Analysis (NB02)
- **Differential abundance:** Mann-Whitney U test, Benjamini-Hochberg correction (FDR < 0.05)
- **Correlation:** vectorised Spearman rank-matrix multiplication; partial correlations with age/BMI/sex residuals
- **Stage-stratified correlations:** Healthy / Early_CRC / Advanced_CRC subgroups
- **Species co-abundance:** upper-triangle Spearman |ρ| ≥ 0.20 matrix
- **Beta diversity:** Aitchison PCoA + PERMANOVA (999 permutations, scikit-bio)
- **Hub species:** ≥ 3 significant metabolite edges

### 4. Machine Learning (NB03)
- **Target selection:** top 15 dysregulated metabolites (early + advanced, sorted by effect size |r|)
- **Models:** XGBoost, LightGBM, ElasticNet, SVR (RBF kernel), Random Forest
- **Validation:** 10-fold stratified cross-validation (stratified on 3-group CRC stage)
- **SHAP:** TreeExplainer (XGBoost/LightGBM/RF), LinearExplainer (ElasticNet), KernelExplainer (SVR, 50 bg / 100 test samples)
- **LODO:** LeaveOneGroupOut across YACHIDA + KONG_EOCRC_2023 + SINHA_VOGTMANN + BORENSTEIN cohorts (when data available)

### 5. MOFA+ / FDA / Trajectory (NB04)
- **Dimensionality reduction:** MOFA+ (`mofapy2`) or TruncatedSVD fallback for negative CLR values
- **Factor–stage association:** Kruskal-Wallis + BH correction
- **Stage separation:** Fisher's Discriminant Analysis (10-fold CV balanced accuracy)
- **Ordinal regression:** `statsmodels` OrderedModel (proportional-odds logit); McFadden pseudo-R²
- **MP arc:** metachronous polyp recovery trajectory with IQR ribbons

### 6. Mediation Network (NB05)
- **Network:** bipartite species–metabolite partial-correlation pairs; hub species/metabolites by degree
- **Causal mediation:** bootstrap ACME (1000-iteration percentile CI + 500-permutation p-values)
- **Multiple testing:** BH correction across all mediation pairs
- **Centrality:** degree + betweenness (NetworkX)

### 7. GutSMASH Benchmarking (NB06)
- **GutSMASH** automated runner: BGC cluster detection per genome (subprocess, 30-min timeout)
- **KEGG mapping:** 85+ KEGG compound IDs → 20+ GutSMASH BGC metabolite classes (expanded from 3 to cover ~43% of SHAP targets)
- **Scope classification:** `gutsmash_relevant` / `primary_metabolism` / `no_kegg_id` / `unmapped`
- **Novel candidates:** SHAP-significant producers without GutSMASH BGC evidence
- **Literature evidence tiers (Tier 0–4):** automated PubMed E-utilities mining; 5 queries per (genus, metabolite) pair; cached to JSON
- **Composite score:** SHAP importance × (1 + literature_tier / 4)
- **Visualisations:** continuous scatter, producer catalogue heatmap, bipartite network, literature evidence heatmap, novel candidate forest plot, Holy Trinity figure

### 8. Mechanistic Flux (NB09)
- **Community metabolic models:** AGORA103 database; `micom.workflows.build()`
- **Growth simulation:** cooperative tradeoff (λ = 0.5, recommended for shotgun metagenomics); linearised L1 norm (GLPK-compatible)
- **Exchange flux analysis:** `strategy="none"` LP-feasible flux distributions; interpret as presence/absence + relative magnitude (not absolute rates)
- **Heatmap atlas:** sample × metabolite (ClusterGrid), stage × metabolite (absolute + Z-score panels), top-30 species × metabolites
- **SHAP cross-reference:** KEGG-level metabolite overlap with NB03 producer candidates; genus-level confirmation
- **Coverage note:** AGORA103 covers ~6,000 named cultured species; MetaPhlAn4/GTDB has ~57,000 genome bins — ~70% of community by abundance is unmodelled

---

## Key Results

### CRC-Associated Metabolites
- Stage-progressive dysregulation identified in early (Stage 0, HS) and advanced (Stage I/II, III/IV) CRC relative to Healthy controls
- Polyamine pathway metabolites (putrescine, cadaverine, spermidine, agmatine, ornithine) among top dysregulated targets
- SCFA (butyrate, propionate, acetate) and secondary bile acids show complementary stage-dependent patterns

### Microbial Producers (Holy Trinity Framework)
- Statistical (SHAP): top-50 species ranked by predictive contribution per metabolite
- Genomic (GutSMASH): BGC-confirmed producers; novel candidates flagged with Tier 0–4 literature evidence
- Mechanistic (MICOM): flux-confirmed excretion across CRC stage progression

### ML Performance
- 5-model benchmark across top 15 dysregulated metabolites
- Best-performing models provide SHAP-interpretable producer rankings
- LODO cross-cohort validation (in progress, pending additional cohort data)

---

## User Data Pipeline

Anyone can run this pipeline on their own paired metagenomics + metabolomics data using Docker.

### Accepted input formats

| File | Accepted formats |
|------|-----------------|
| **Species table** | MetaPhlAn3/4 profile (clade names with `k__\|s__`), generic OTU/ASV TSV, or any tab/comma-separated abundance matrix (auto-transposed if needed) |
| **Metabolomics table** | KEGG IDs (`C12345`), HMDB IDs (`HMDB00062`), common compound names (putrescine, leucine, …), or generic feature names |
| **Metadata** | Any TSV/CSV with a `sample_id` / `Sample` / `SampleID` column plus at least one grouping/condition column |

### Running with Docker

```bash
# 1. Clone the repository
git clone https://github.com/axn14/Correlating-the-MGx-MBx-Axis.git
cd Correlating-the-MGx-MBx-Axis

# 2. Build and start (first build takes ~10 min to install all dependencies)
docker-compose up --build

# 3. Open http://localhost:8501 in your browser
# 4. Click "Upload & Run" in the sidebar and follow the on-screen steps
```

Results are written to `./results/` on your host machine and are immediately available in the viewer pages.

### What runs on user data

| Step | Available | Notes |
|------|-----------|-------|
| Preprocessing (NB01) | ✅ | CLR + log₁₀ transforms, prevalence filters |
| Statistical associations (NB02) | ✅ | DA, Spearman, partial corr (auto-detects confounders) |
| ML benchmarking (NB03) | ✅ | 5-model nested CV; 10 Optuna trials (reduced for speed) |
| MOFA+ (NB04) | ✅ | Requires R in container; gracefully skipped if unavailable |
| Network & mediation (NB05) | ✅ | Partial-corr network + bootstrap ACME |
| GutSMASH BGC (NB06, E8) | ❌ | Requires assembled MAGs — not generalisable to abundance tables |
| MICOM flux (NB09, E9) | ❌ | Requires AGORA103 model coverage for user species |

**Evidence streams for user data:** E1–E7 only. TRINITY criterion for user data = E1/E2 (statistical) + E3 (KEGG enzyme) + E7 (within-stage Spearman).

---

## Limitations

- **Not clinically validated:** This tool has been developed and evaluated on the YACHIDA-CRC-2019 cohort only. Do not use results for clinical or diagnostic decisions.
- **E8/E9 not available for user data:** GutSMASH (E8) requires assembled metagenome-assembled genomes (MAGs); MICOM (E9) requires your species to be present in the AGORA103 database. Both are therefore unavailable when submitting abundance tables.
- **ML runtime:** The ML step (step 3) runs up to 45 targets × 10 outer folds × 5 models and may take 20–60 minutes depending on dataset size and CPU.
- **MOFA+ requires R:** The Dockerfile installs R base. If `mofapy2` fails to initialise, MOFA+ is automatically skipped with a warning.
- **AGORA103 coverage:** ~70% of MetaPhlAn4 species bins are not in AGORA103 (which covers ~6,000 cultured species). MICOM results for user data would reflect only the covered fraction.

---

## Installation

### Recommended: conda environment
```bash
conda env create -f environment.yml
conda activate metabolomics-pipeline
```

### pip (Python 3.12)
```bash
pip install -r requirements.txt
```

### Additional tools (upstream pipeline, NB00 only — requires Linux/WSL)
```bash
# MetaPhlAn4 + HUMAnN3 + fastp via bioconda (see environment.yml)
conda install -c bioconda metaphlan humann fastp
```

### MICOM (NB09)
```bash
pip install micom
# Download AGORA103 database:
# https://zenodo.org/record/7548056 → agora103_species.qza → Data/micom/
# Western diet medium:
# https://github.com/micom-dev/media → western_diet_gut.csv → Data/micom/
```

### GutSMASH (NB06)
```bash
pip install gutsmash          # or see appendix cell in NB06 for conda install
# GutSMASH genomes → Data/genomes/{genome}.fna
```

---

## Data

Raw data files are not tracked in this repository. Place files as:
```
Data/
├── {DATASET} metadata.tsv
├── {DATASET} species.tsv
├── {DATASET} mtb.tsv
├── {DATASET} mtb.map.tsv
├── micom/
│   ├── agora103_species.qza
│   └── western_diet_gut.csv
└── genomes/
    └── *.fna   (for GutSMASH, NB06)
```

Datasets used: `YACHIDA-CRC-2019`, `SINHA_CRC_2016`

Results directory (figures, tables, pickles) is also excluded from tracking.

---

## Dependencies

See `requirements.txt` (pip) and `environment.yml` (conda) for full pinned dependency lists.

Core packages: `pandas`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`,
`statsmodels`, `networkx`, `xgboost`, `lightgbm`, `shap`, `micom`, `mofapy2`,
`pyopenms`, `requests` (PubMed API), `scikit-bio` (PERMANOVA)

---

## Citation

If you use this pipeline, please cite:

- **MICOM:** Diener et al. (2020). MICOM: Metagenome-Scale Modeling To Infer Metabolic Interactions in the Gut Microbiota. *mSystems*, 5(1).
- **GutSMASH:** Martínez-Martínez et al. (2024). GutSMASH predicts specialized primary metabolic gene clusters in the gut microbiota. *Nature Methods*.
- **AGORA103:** Heinken et al. (2023). Genome-scale metabolic reconstruction of 7,302 human microorganisms for personalized medicine. *Nature Biotechnology*.
- **Yachida et al. (2019).** Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the gut microbiota in colorectal cancer. *Nature Medicine*, 25, 968–976.

---

## License

For academic use. Contact repository owner for commercial use enquiries.
