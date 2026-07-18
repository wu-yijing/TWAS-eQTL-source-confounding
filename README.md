# TWAS eQTL Source Confounding — Systematic Evaluation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21428347.svg)](https://doi.org/10.5281/zenodo.21428347)

## Overview

This repository contains analysis scripts and processed data for:

**"eQTL source confounding in TWAS: a dual-source 2×2 decomposition"**

Target journal: **Genetic Epidemiology** (under consideration)

> The original v1.0.0 release (Zenodo [10.5281/zenodo.21428347](https://doi.org/10.5281/zenodo.21428347)) corresponded to the iScience submission titled *"eQTL source confounding systematically biases TWAS cross-population replication …"* (HOTAIR binding proteins, 104 genes, three diabetic complications). v2.0.0 adds the schizophrenia (SCZ) dual-source 2×2 decomposition replication and retargets the manuscript to Genetic Epidemiology.

### Core Question

Does the choice of eQTL weight source (GTEx v8 tissue-specific vs. eQTLGen large-scale whole blood) systematically alter TWAS cross-population replication conclusions?

### Key Findings

| Metric | Value |
|--------|-------|
| Total genes analyzed | 104 (30 candidate + 44 non-candidate + 30 T2DM control) |
| Spearman ρ (GTEx vs eQTLGen Z-scores) | 0.289 |
| Direction consistency (all pairs) | 59.5% |
| Direction consistency (candidate genes only) | 37.5% |
| GTEx FDR rate (candidate) | 40.7% |
| eQTLGen FDR rate (candidate) | 0% |
| RNH1 GTEx Z-score | +8.50 (FDR q=1.5×10⁻¹⁵) |
| RNH1 eQTLGen Z-score | +1.18 (P=0.237) |
| Cross-population I² (RNH1, 3 cohorts) | 95.0% |

## SCZ 2×2 Decomposition Replication (added in v2.0.0)

A dual-source 2×2 decomposition was applied to schizophrenia (SCZ, PGC3 wave3 European) to test whether the source/panel vs tissue-context structure generalizes beyond the metabolic discovery cohort.

- **Method**: self-implemented TWAS Z = (w′R⁻¹z)/√(w′R⁻¹w) with LD from 1000G EUR, ridge RIDGE=0.1; 2×2 axes = source (eQTLGen Whole_Blood vs GTEx multi-tissue Stouffer) and tissue (GTEx Whole_Blood vs GTEx Nerve_Tibial).
- **Results (n = 2,511 gene–tissue pairs)**:
  - Source axis ρ = +0.522 (68.3% direction concordance)
  - Tissue axis ρ = +0.509 (71.2% direction concordance)
  - Both axes positive and comparable in magnitude — the axis ordering is **trait-dependent**, not a fixed "source always dominates" rule.
- **Sparsity robustness**: GTEx v8 MASHR models are sparse (median 2 SNPs/gene). Stratifying by model size (2 / 3–4 / ≥5 SNPs) gave stable ρ (≈0.5 across strata) with bootstrap 95% CI excluding zero; the eQTLGen side (median 786 SNPs/gene) anchors the source axis. See `scz_replication/results/`.
- **Scripts**: `scz_replication/scz_twas.py` (TWAS + decomposition), `scz_replication/build_weights_db.py` (weight DB), `scz_replication/robustness_scz.py` (stratified robustness), `scz_replication/make_figure9_scz.py` + `scz_replication/make_robustness_fig.py` (figures).
- **Figures**: `figures/scz/Fig9_2x2_decomposition.*` (source/tissue scatter) and `figures/scz/FigS1_sparsity_robustness.*` (stratified robustness).

## Quick Start (Docker — Recommended)

```bash
# 1. Build the Docker image
docker build -t twas-eqtl-repro .

# 2. Run the full reproducibility pipeline
docker run --rm -v $(pwd)/output:/app/output twas-eqtl-repro

# 3. View outputs
open output/figs/    # All figures in PDF + PNG
open output/logs/    # Analysis provenance logs
```

### What the Docker image contains
- **Python 3.13 + R 4.5.2** with all dependencies (pandas, scipy, statsmodels, MatchIt, etc.)
- **MetaXcan** (for S-PrediXcan, run separately with `run_spredixcan.sh`)
- **Downstream analysis scripts** that generate all paper figures
- **Processed data** needed to reproduce figures from intermediate results

## Local Setup (Conda)

```bash
# With your local conda installation:
conda env create -f environment.yml
conda activate twas-eqtl

# Run downstream analysis
python scripts/python/04_generate_all_figures.py
```

## Reproducing the Full Pipeline (Optional)

To reproduce from raw GWAS data (requires large input files):

```bash
# Mount external data and run S-PrediXcan
docker run --rm \
  -v /path/to/finngen_data:/input/finngen \
  -v /path/to/gtex_models:/input/gtex \
  -v /path/to/eqtlgen_weights:/input/eqtlgen \
  -v /path/to/1000g_ref:/input/1000g \
  -v $(pwd)/output:/app/output \
  twas-eqtl-repro \
  bash /app/run_spredixcan.sh
```

## Repository Structure

```
├── Dockerfile                     # Container definition
├── environment.yml                # Conda environment (pinned)
├── requirements.txt               # pip dependencies
├── renv.lock                      # R environment lock
├── run_all.sh                     # Main entry point (Docker CMD)
├── run_spredixcan.sh              # S-PrediXcan execution reference
├── .dockerignore                  # Docker build exclusions
├── data/
│   ├── raw/                       # (empty — raw data too large for repo)
│   └── processed/                 # Intermediate results (CSV)
├── scripts/
│   ├── python/                    # Python analysis scripts
│   └── R/                         # R scripts (Mahalanobis matching)
├── figs/                          # Generated figures
├── docs/                          # Documentation
└── analyses/logs/                 # Provenance logs + SHA256 checksums
```

```
TWAS-eQTL-source-confounding/
├── scripts/
│   ├── python/
│   │   ├── 01_visualize_mahalanobis_love_plot.py    # Love plot (covariate balance)
│   │   ├── 02_density_scatter_consistency.py         # Density, scatter, bar charts
│   │   ├── 03_enrichment_analysis.py                 # Enrichment + statistical tests
│   │   ├── 04_generate_all_figures.py                # All main + supp figures
│   │   ├── 05_generate_decision_framework.py         # Decision framework flowchart (Fig 6)
│   │   ├── 06_generate_supplementary_figures.py      # Supplementary figures S5-S7
│   │   ├── 07_generate_supplementary_tables.py       # Supplementary tables S4-S6
│   │   └── 08_generate_tables_S1_S2.py               # Supplementary tables S1-S2
│   └── R/
│       └── (placeholder)                             # matchit R scripts
├── data/
│   ├── raw/                                          # (public data references only)
│   └── processed/
│       ├── eqtlgen_vs_gtex_comparison.csv            # GTEx vs eQTLGen Z-score pairs
│       ├── eqtlgen_spredixcan_results.csv             # eQTLGen TWAS full results
│       ├── mahalanobis_matched_pairs.csv              # Mahalanobis matching pairs
│       ├── matched_enrichment_results.csv             # Post-matching enrichment
│       ├── enrichment_comparison.csv                  # Enrichment by group × phenotype × eQTL
│       ├── covariate_matrix.csv                       # Gene-level covariates
│       ├── layer_analysis.csv                         # Pull-down vs literature stratification
│       ├── viz_z_distribution.csv                     # Z distribution data for viz
│       └── candidate_comparison_DR.csv                # Candidate DR comparison data
├── figs/                                              # Output figures (generated)
├── docs/                          # Documentation
├── environment.yml                                    # Conda environment
├── requirements.txt                                   # Python dependencies
├── .gitignore
├── LICENSE (MIT)
└── README.md
```

## SCZ Replication Folder Layout

```
scz_replication/
├── scz_twas.py                  # TWAS Z + 2×2 decomposition (self-implemented)
├── build_weights_db.py          # assemble eqtlgen / GTEx WB / GTEx NT weight DB
├── robustness_scz.py           # sparsity-stratified robustness + bootstrap CI
├── make_figure9_scz.py          # Fig9 source/tissue scatter
├── make_robustness_fig.py       # FigS1 stratified bar chart
├── patch_manuscript_scz.py      # manuscript write-back (SCZ section)
├── patch_robustness.py          # manuscript write-back (limitation)
└── results/
    ├── scz_decomp_limit0.json  # decomposition ρ / n / concordance
    ├── scz_robustness.json     # stratified robustness + CI
    └── scz_twas_results_limit0.csv  # per-gene TWAS Z (10,357 genes)
```
Manuscript (GE-retitled) and cover letter: `manuscript/`. SCZ figures: `figures/scz/`.

> Large intermediates (`weights.db`, extracted eQTLGen RDat weights, raw GWAS) are excluded by `.gitignore`; regenerate via the scripts above.

## Data Sources

All GWAS and eQTL summary statistics used in this study are from **publicly available sources**:

| Dataset | Source | Access |
|---------|--------|--------|
| FinnGen R13 (DR/DN/DPN) | [FinnGen](https://www.finngen.fi/) | Release 13 |
| GTEx v8 MASHR weights | [GTEx Portal](https://gtexportal.org/) | Public access |
| eQTLGen whole blood cis-eQTL | [eQTLGen](https://www.eqtlgen.org/) | Public access |
| 1000 Genomes EUR LD | [1000 Genomes](https://www.internationalgenome.org/) | Public access |
| UK Biobank DR (GCST90043640) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) | Public access |
| UK Biobank DR (Xue et al. 2022) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) | Public access |
| Cai et al. 2026 UKB T2D-DR | (in preparation) | Apply to UKB |
| Sakaue et al. 2021 DN (ebi-a-GCST90018832) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) | Public access |

## Requirements

### Python
```bash
pip install -r requirements.txt
```

Or with conda:
```bash
conda env create -f environment.yml
conda activate twas-eqtl
```

### R
- `matchit` (>= 4.5.0) — for Mahalanobis matching
- S-PrediXcan (v1.0, MetaXcan v0.8.1) — see [https://github.com/hakyimlab/MetaXcan](https://github.com/hakyimlab/MetaXcan)

## Usage

```bash
# Step 1: Run enrichment analysis and statistical tests
python scripts/python/03_enrichment_analysis.py

# Step 2: Generate figures
python scripts/python/01_visualize_mahalanobis_love_plot.py
python scripts/python/02_density_scatter_consistency.py
```

## Reproducibility

All processed data tables are provided in `data/processed/` and `scz_replication/results/`. Analysis scripts are version-controlled in this repository under MIT license. The repository snapshot is archived at Zenodo: [https://doi.org/10.5281/zenodo.21428347](https://doi.org/10.5281/zenodo.21428347). A new Zenodo version (v2.0.0, including the SCZ replication) is minted automatically when the corresponding GitHub Release is published.

## Citation

Wu Y, Chen M, Wu Q, Zhao J, Jin G. (2026). TWAS eQTL Source Confounding — Systematic Evaluation (v1.0.0). Zenodo. [https://doi.org/10.5281/zenodo.21428347](https://doi.org/10.5281/zenodo.21428347)

Wu Y, Chen M, Wu Q, Zhao J, Jin G. (2026). TWAS eQTL Source Confounding — Systematic Evaluation (v2.0.0, SCZ 2×2 decomposition replication, retargeted to Genetic Epidemiology). Zenodo. [https://doi.org/10.5281/zenodo.21428347](https://doi.org/10.5281/zenodo.21428347)

## License

MIT License. See [LICENSE](LICENSE) for details.
