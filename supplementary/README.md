# Supplementary Information — README

**Manuscript:** *eQTL source selection systematically biases TWAS cross-population replication — quantitative evidence from HOTAIR binding proteins across 104 genes and three diabetic complications.*

**Target journal:** iScience (Cell Press)

This folder contains the supplementary tables and figures referenced in the main
manuscript (`manuscript_iScience_v2.docx`). The main-text novelty figure
**Figure 8 (2×2 discordance decomposition)** is part of the main manuscript and is
*not* duplicated here.

---

## 1. File-naming convention

Every supplementary item is provided as a triple:

| Suffix | Format | Purpose |
|--------|--------|---------|
| `.pdf` | vector (publication-quality) | final submission / printing |
| `.png` | raster preview (300 dpi) | quick review / inline |
| `.txt`  | plain-text legend | caption text (also embedded in manuscript) |

Tabular items add data files: `Table_1`–`Table_3` ship as `.csv`;
`Table_S1`–`Table_S5` ship as both `.csv` and `.xlsx`.

---

## 2. Main supplementary tables (cited as **Table 1–3** in the main text)

| File stem | Content |
|-----------|---------|
| `Table_1_gtex_vs_eqtlgen_comparison` | GTEx vs eQTLGen Z-score comparison for core genes. Extreme eQTL-source-dependent discordance across DR/DN/DPN (e.g. RNH1 in DR: GTEx Z = +13.82 vs eQTLGen Z = +1.18). |
| `Table_2_fdr_enrichment_rates` | FDR-significant enrichment rates (%) by gene group (Candidate / Non-Candidate / T2DM Control), phenotype (DR/DN/DPN), and eQTL source (GTEx v8 / eQTLGen). Shows candidate enrichment inversion (GTEx 40.7% → eQTLGen 0%). |
| `Table_3_integrated_evidence_assessment` | Integrated evidence assessment across all analytical layers for the primary claims (RNH1 top candidate; candidate enrichment; cross-population replicability) with overall verdicts. |

## 3. Supplementary tables (**Table S1–S5**)

| File stem | Content |
|-----------|---------|
| `Table_S1_Supplementary_Material` | Gene list with group assignment (Candidate / NonCandidate / T2DM_Control), identification source, RNA pull-down MS metrics (Unused score, sequence coverage, peptides), and genomic covariates (length, GC%, mean eQTL SNP count). |
| `Table_S2_Supplementary_Material` | GTEx v8 baseline TWAS results for all analyzed genes across Nerve_Tibial & Whole_Blood and DR/DN/DPN: Z, P, model SNP counts, Stouffer multi-tissue Z (`Z_Multi`), ACAT-O (`P_ACAT_O`), and FDR q-values. |
| `Table_S3_Supplementary_Material` | Gene covariate matrix for all 104+ genes (30 candidates, 44+ non-candidates, 30 T2DM controls): log10 length, Length_bp, GC%, eQTL SNP count, pull-down Unused score. |
| `Table_S4_Supplementary_Material` | Mahalanobis matched pairs — 30 nearest-neighbour pairs (one candidate + one matched non-candidate) with covariates used for matching (log10 length, GC%, eQTL SNP count). |
| `Table_S5_Supplementary_Material` | Cross-population replication results for RNH1 and key candidates across FinnGen R13, UKB (GCST90043640; Xue et al. 2022), and Sakaue et al. 2021 (DN), with meta-analysis and heterogeneity (Cochran Q = 40.0, I² = 95.0%, 95% PI [−4.00, +10.68]). |

## 4. Supplementary figures (**Figure S1–S8**)

| File stem | Content |
|-----------|---------|
| `Figure_S1_Supplementary_Material` | Mahalanobis matched pairs — covariate distribution (log10 length, eQTL SNP count, GC%) with candidate–control linkage lines; confirms balanced matching. |
| `Figure_S2_Supplementary_Material` | Pull-down vs literature stratification of eQTLGen FDR rates across DR/DN/DPN by identification source. |
| `Figure_S3_Supplementary_Material` | Multi-tissue integration method sensitivity: Stouffer Z-meta (15 pairs) vs ACAT-O (14 pairs) — high concordance. |
| `Figure_S4_Supplementary_Material` | \|Z\| density by gene group (eQTLGen whole blood) across DR/DN/DPN; \|Z\| = 1.96 reference. |
| `Figure_S5_Supplementary_Material` | eQTL SNP count distribution: per-group boxplot (Panel A) and GTEx v8 vs eQTLGen cis-eQTL density (~20-fold difference, Panel B). |
| `Figure_S6_Supplementary_Material` | Direction-consistency sensitivity across \|Z\| thresholds (0.01–0.50); stable ~60% rate vs 50% random expectation. |
| `Figure_S7_Supplementary_Material` | \|Z\| distribution by gene group (eQTLGen, log10 scale); non-candidates show strongest signals. |
| `Figure_S8_Supplementary_Material` | Enrichment comparison bar chart of FDR rates (DR/DN/DPN × three groups) under GTEx v8 (red) vs eQTLGen (blue). |

---

## 5. Data availability

Processed data tables and analysis scripts are archived at **Zenodo**:
<https://doi.org/10.5281/zenodo.21238203>

The live, version-controlled code repository (this repository) is released under
the MIT license. Containerized reproducibility is provided via the project
`Dockerfile` (see repository root).

---

## 6. Pre-submission checklist (author action required)

- [x] **Legend numbering drift — FIXED.** The `Table_S3` / `Table_S4` /
      `Table_S5` legend headers previously read "Supplementary Table S4" /
      "S5" / "S6"; they are now corrected to **S3 / S4 / S5** to match the file
      names and the main-text citations.
- [x] **Figure S8 legend header — FIXED.** `Figure_S8_..._txt` now begins with the
      "Supplementary Figure S8. Enrichment Comparison Across Gene Groups and eQTL
      Sources." header line, consistent with S1–S7.
- [ ] Confirm every `Figure S#` legend text matches the corresponding image
      content (cross-check against the manuscript figure-legend list).
- [ ] Verify `Table S1–S5` / `Figure S1–S8` numbering in the manuscript matches
      these file names exactly.
