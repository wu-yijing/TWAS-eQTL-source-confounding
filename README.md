# TWAS eQTL Source Confounding — Systematic Evaluation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21238202.svg)](https://doi.org/10.5281/zenodo.21238202)

## Overview

This repository contains analysis scripts and processed data for:

**"eQTL-source discordance in TWAS: a dual-source sensitivity analysis with true-negative calibration"**

Status: code & data release accompanying the manuscript (BMC Genomics submission).

> The original v1.0.0 release (Zenodo [10.5281/zenodo.21428347](https://doi.org/10.5281/zenodo.21428347)) corresponded to the iScience submission titled *"eQTL source confounding systematically biases TWAS cross-population replication …"* (HOTAIR binding proteins, 104 genes, three diabetic complications). The archive now resolves through the concept DOI [10.5281/zenodo.21238202](https://doi.org/10.5281/zenodo.21238202), which always points to the latest version; the current version is v2.4.0 ([10.5281/zenodo.22691102](https://doi.org/10.5281/zenodo.22691102)), adding the gene-level cluster-robustness re-analysis (S1) on top of the SCZ dual-source decomposition replication.

### Core Question

Does the choice of eQTL weight source (GTEx v8 tissue-specific vs. eQTLGen large-scale whole blood) systematically alter TWAS cross-population replication conclusions?

### Key Findings

> **Note (2026-09-09, harmonized recompute):** The table below reflects the **current manuscript version (BMC Genomics submission)**, in which all eQTLGen-arm values were recomputed after **three-way allele harmonization** (eQTLGen assessed/other alleles ↔ FinnGen R13 ref/alt ↔ 1000G EUR A2) using the official eQTLGen FDR < 0.05 weight catalogue. This supersedes both the v1.0.0 (iScience) pipeline values and the earlier archived eQTLGen export, which contained an allele-alignment defect (24.1% of Z-scores flipped sign after harmonization). FDR enrichment rates are computed by Benjamini–Hochberg correction within each group × phenotype × source stratum; 41 of 144 genes lacked an eQTLGen model in the official FDR < 0.05 catalogue and are disclosed in the manuscript Methods. See `data/processed/eqtlgen_spredixcan_harmonized_results.csv` (gene-level Z-scores) and `data/processed/enrichment_comparison_harmonized.csv` (stratum-wise enrichment rates with Clopper–Pearson CIs).

| Metric | Value | Notes |
|--------|-------|-------|
| Total genes analyzed | 104 (30 candidate + 44 non-candidate + 30 T2DM control) | 114 records in covariate matrix (includes variants) |
| Spearman ρ (GTEx WB vs harmonized eQTLGen Z-scores) | 0.32 (P = 0.001) | 96 gene–phenotype pairs from 32 genes (primary arm) |
| Direction consistency (primary arm) | 63.5% (61/96) | Exact binomial P = 0.010 vs 50%; DR 65.6%, DN 65.6%, DPN 59.4% |
| Direction consistency (candidate subset, GTEx NT vs eQTLGen) | 70.0% (42/60) | 20 genes; subset definitions shift estimates by ≤7 pp (manuscript Table 2c) |
| GTEx FDR rate (candidate, DR) | 48.1% (13/27) | GTEx v8 MASHR, multi-tissue Stouffer integration |
| eQTLGen FDR rate (candidate, DR) | 41.7% (10/24) | Harmonized recompute; group rates 41.7–66.7% across phenotypes, comparable to non-candidate and T2DM controls |
| RNH1 GTEx Z-score (DR, Nerve_Tibial) | +13.82 | S-PrediXcan; from sparse 3-SNP MASHR model |
| RNH1 eQTLGen Z-score (DR) | +13.32 | Harmonized recompute; source-stable signal |
| Cross-population heterogeneity (RNH1, DR; FinnGen R13 vs UK Biobank) | I² = 98.7% | Random-effects pooled Z = +7.13, 95% prediction interval −13.77 to +28.04 (spans zero) |
| SCZ replication (n=2,511 pairs) | 68.3% direction consistency | Axis difference not significant (bootstrap P = 0.530) |

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
│   │   ├── 08_generate_tables_S1_S2.py               # Supplementary tables S1-S2
│   │   └── s1_cluster_robustness/                    # Gene-level cluster-robust re-analysis (S1)
│   │       ├── s1_recon.py                           # Reconstruct analysis arms from raw tables
│   │       ├── s1_diag.py                            # Diagnose NaN genes / cluster sizes
│   │       ├── s1_fix.py                             # Confirm NaN genes, rebuild 96-pair arm
│   │       ├── s1_cluster.py                         # ICC/DEFF, cluster-robust variance, bootstraps
│   │       ├── s1_scz_cluster.py                     # Genome-wide SCZ arm re-analysis
│   │       ├── s1_results.json                       # HOTAIR arms (machine-readable)
│   │       ├── s1_scz_results.json                   # SCZ arm (machine-readable)
│   │       ├── s1_primary_arm_96pairs.csv            # 96 primary-arm pairs (Z + direction)
│   │       └── s1_anchor_102pairs.csv                # 102 anchor-set pairs
│   └── R/
│       └── (placeholder)                             # matchit R scripts
├── data/
│   ├── raw/                                          # (public data references only)
│   └── processed/
│       ├── eqtlgen_vs_gtex_comparison.csv            # GTEx vs eQTLGen Z-score pairs
│       │                                             #   ⚠️ DEPRECATED for direction-consistency analyses:
│       │                                             #   predates the three-way allele harmonisation, so its
│       │                                             #   Same_Direction column reflects pre-harmonisation Z.
│       │                                             #   Use eqtlgen_spredixcan_harmonized_results.csv instead.
│       ├── eqtlgen_spredixcan_results.csv             # eQTLGen TWAS results (archived export; superseded — see harmonized file below)
│       ├── eqtlgen_spredixcan_harmonized_results.csv  # eQTLGen TWAS results, three-way allele-harmonized recompute (current manuscript)
│       ├── enrichment_comparison_harmonized.csv       # Stratum-wise enrichment rates (harmonized), Clopper–Pearson CIs
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
Manuscript and cover letter: `manuscript/`. SCZ figures: `figures/scz/`.

> Large intermediates (`weights.db`, extracted eQTLGen RDat weights, raw GWAS) are excluded by `.gitignore`; regenerate via the scripts above.

## Cluster-Robustness Re-analysis of Direction Consistency (S1)

The primary between-source comparison comprises 96 gene–phenotype pairs that derive from
only **32 genes**, so a gene's three phenotypes (DR, DN, DPN) are not independent
observations. The `scripts/python/s1_cluster_robustness/` folder quantifies this clustering
and recomputes the direction-consistency test with three cluster-aware procedures.

**Result: the intra-gene correlation is negligible (ICC = 0.000, design effect = 1.000).**
The exact binomial test reported in the manuscript (P = 0.010) is therefore **not materially
anti-conservative** — a cluster-robust sandwich test gives P = 0.013, a gene-level cluster
bootstrap gives P = 0.017 (95% CI 53.1–74.0%), and a gene-label permutation test gives
P = 0.008. All variants remain significant at α = 0.05, and the CI widens by only 0.9
percentage points (53.1–74.0% vs 53.1–73.1%).

For the genome-wide SCZ arm the complete-case set contributes **exactly one pair per gene**
(2,511 pairs from 2,511 genes), so the design effect is 1 by construction and no cluster
correction is required; the gene-resampling bootstrap coincides with the exact binomial test
(68.26%, ρ = 0.522 for the resource/sample-size axis; 71.21%, ρ = 0.509 for the
tissue-context axis).

Seeds: `20260910` (HOTAIR arms) and `20260726` (SCZ arm), B = 10,000 each.
Full methods, reproduction steps and a stale-file caveat are in
`scripts/python/s1_cluster_robustness/README.md`.

## Data Sources

All GWAS and eQTL summary statistics used in this study are from **publicly available sources**:

| Dataset | Source | Access |
|---------|--------|--------|
| FinnGen R13 (DR/DN/DPN) | [FinnGen](https://www.finngen.fi/) | Release 13 |
| GTEx v8 MASHR weights | [GTEx Portal](https://gtexportal.org/) | Public access |
| eQTLGen whole blood cis-eQTL | [eQTLGen](https://www.eqtlgen.org/) | Public access |
| 1000 Genomes EUR LD | [1000 Genomes](https://www.internationalgenome.org/) | Public access |
| UK Biobank DR (GCST90043640) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) | Public access |
| UK Biobank DR (Xue et al. 2022; ieu-b-4803, 14,147 cases/322,390 controls) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) | Public access |
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

All processed data tables are provided in `data/processed/` and `scz_replication/results/`. Analysis scripts are version-controlled in this repository under MIT license. The repository snapshot is archived at Zenodo under the concept DOI [https://doi.org/10.5281/zenodo.21238202](https://doi.org/10.5281/zenodo.21238202), which always resolves to the latest version (current version v2.4.0, [10.5281/zenodo.22691102](https://doi.org/10.5281/zenodo.22691102)). The gene-level cluster-robustness re-analysis (S1) is archived both here (`scripts/python/s1_cluster_robustness/`) and in that Zenodo release as `s1_cluster_robustness.zip`.

> ⚠️ `data/processed/eqtlgen_vs_gtex_comparison.csv` is **superseded for direction-consistency analyses**: it predates the three-way allele harmonisation, so its `Same_Direction` column reflects pre-harmonisation Z-scores. Use `data/processed/eqtlgen_spredixcan_harmonized_results.csv` instead.

## Citation

Wu Y, Chen M, Wu Q, Zhao J, Jin G. (2026). TWAS eQTL Source Confounding — Systematic Evaluation (v1.0.0). Zenodo. [https://doi.org/10.5281/zenodo.21428347](https://doi.org/10.5281/zenodo.21428347)

Wu Y. (2026). eQTL-source discordance in TWAS: a dual-source sensitivity analysis with true-negative calibration (code & data) (v2.4.0). Zenodo. [https://doi.org/10.5281/zenodo.22691102](https://doi.org/10.5281/zenodo.22691102)

## License

MIT License. See [LICENSE](LICENSE) for details.
