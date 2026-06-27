# TWAS eQTL Source Confounding — Systematic Evaluation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository contains analysis scripts and processed data for:

**"eQTL source confounding systematically biases TWAS cross-population replication: a quantitative evaluation in HOTAIR binding proteins across 104 genes and three diabetic complications"**

Target journal: **GigaScience**

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

## Repository Structure

```
TWAS-eQTL-source-confounding/
├── scripts/
│   ├── python/
│   │   ├── 01_visualize_mahalanobis_love_plot.py    # Love plot (covariate balance)
│   │   ├── 02_density_scatter_consistency.py         # Density, scatter, bar charts
│   │   └── 03_enrichment_analysis.py                 # Enrichment + statistical tests
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
├── docs/
│   └── gigascience_paper_outline.md                   # Restructured paper outline
├── environment.yml                                    # Conda environment
├── requirements.txt                                   # Python dependencies
├── .gitignore
├── LICENSE (MIT)
└── README.md
```

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
- S-PrediXcan (v1.0) — see [https://github.com/hakyimlab/MetaXcan](https://github.com/hakyimlab/MetaXcan)

## Usage

```bash
# Step 1: Run enrichment analysis and statistical tests
python scripts/python/03_enrichment_analysis.py

# Step 2: Generate figures
python scripts/python/01_visualize_mahalanobis_love_plot.py
python scripts/python/02_density_scatter_consistency.py
```

## Reproducibility

All processed data tables are deposited in **GigaDB** under accession [DOI pending]. Analysis scripts are version-controlled in this repository under MIT license.

## Citation

[Citation information — to be added upon publication]

## License

MIT License. See [LICENSE](LICENSE) for details.
