# Figure Alt Text — GigaScience Submission
# 提交前请仔细核对每张图的 Alt Text 描述是否准确

## Figure 1 — GTEx vs eQTLGen TWAS Z-score Comparison
**File:** Fig1_Scatter_Plot.pdf (166 × 60mm)
**Alt text:** Scatter plot comparing TWAS Z-scores from GTEx v8 and eQTLGen across three diabetic complications (DR, DN, DPN), showing low concordance between eQTL sources. Circles indicate same-direction genes; crosses indicate opposite-direction genes. Spearman ρ = 0.289 (P < 0.001) overall.

## Figure 2 — Core Gene |Z-score|: GTEx vs eQTLGen (DR)
**File:** Fig2_Waterfall.pdf (166 × 82mm)
**Alt text:** Horizontal bar plot comparing absolute TWAS Z-scores for 11 core genes in DR under GTEx v8 (light red) and eQTLGen (light blue). RNH1 shows the most extreme attenuation from |Z|=8.50 (GTEx) to |Z|=1.18 (eQTLGen). Dashed red line at |Z|=1.96 marks nominal significance threshold.

## Figure 3 — FDR Enrichment Rate: GTEx v8 vs eQTLGen (DR)
**File:** Fig3_Enrichment.pdf (166 × 96mm)
**Alt text:** Grouped bar chart comparing FDR enrichment rates across three gene groups (Candidate, Non-Candidate, T2DM Control) under GTEx v8 (beige) and eQTLGen (light blue). Candidate enrichment drops from 40.7% under GTEx to 0% under eQTLGen.

## Figure 4 — Covariate Balance (Love Plot)
**File:** Fig4_Love_Plot.pdf (94 × 50mm, half-page)
**Alt text:** Love plot showing standardized mean differences for three covariates (gene length, GC content, eQTL SNP count) before (red circles) and after (blue squares) Mahalanobis distance matching. All |SMD| values are below 0.25 after matching.

## Figure 5 — RNH1 Effect Attenuation Across eQTL Sources and Populations
**File:** Fig5_RNH1_Attenuation.pdf (181 × 69mm)
**Alt text:** Two-panel figure. Panel (a): Grouped bar chart showing RNH1 |Z|-scores across DR, DN, and DPN under GTEx v8 (light red) and eQTLGen (light blue), demonstrating dramatic attenuation under the large-sample source. Panel (b): Forest plot showing RNH1 cross-population Z-scores across FinnGen and UK Biobank cohorts with 95% confidence intervals, demonstrating extreme heterogeneity (I² = 95%).

## Figure 6 — TWAS eQTL Source Selection Decision Framework
**File:** Fig6_Decision_Framework.pdf (353 × 199mm, to be scaled to 170mm)
**Alt text:** Flowchart guiding interpretation of TWAS findings based on availability of multiple independent eQTL sources. If multiple sources are available, compute direction consistency and Spearman ρ. Thresholds: ρ > 80% = robust; ρ = 60-80% = caveat needed; ρ < 60% = source-dependent. For single-source studies, mandatory caveats about tissue/sample-size limitations apply. Post-hoc validation includes Mahalanobis matching, cross-population replication, and functional annotation stratification.

## Figure 7 — Study Design Flowchart
**File:** Fig7_Flowchart.pdf (191 × 138mm)
**Alt text:** Schematic overview of the six-module analytical pipeline. Steps include: RNA pull-down MS (69 proteins), public eQTL/GWAS data sources (FinnGen R13, GTEx v8, eQTLGen), gene group assignment (30 candidates, 44 non-candidates, 30 T2DM controls), S-PrediXcan TWAS analysis under both GTEx and eQTLGen, systematic comparison (direction consistency, Spearman ρ, enrichment rate, Mahalanobis matching), and conclusion that eQTL source selection systematically biases TWAS replication.

---

# GigaScience 提交核对清单

## 图片文件
- [ ] Fig1_Scatter_Plot.pdf — 主图 1，矢量格式 ✓
- [ ] Fig2_Waterfall.pdf — 主图 2，矢量格式 ✓
- [ ] Fig3_Enrichment.pdf — 主图 3，矢量格式 ✓
- [ ] Fig4_Love_Plot.pdf — 主图 4，矢量格式（半页宽）✓
- [ ] Fig5_RNH1_Attenuation.pdf — 主图 5，矢量格式 ✓
- [ ] Fig6_Decision_Framework.pdf — 主图 6，矢量格式 ✓
- [ ] Fig7_Flowchart.pdf — 主图 7，矢量格式 ✓
- [ ] 所有 .png 后备格式（非强制，建议保留）✓

## GigaScience 格式要求逐项核对
- [ ] 矢量格式（PDF/eps/SVG）— 全部为 PDF ✓
- [ ] 分辨率 300 DPI — 脚本设置 savefig.dpi=300 ✓
- [ ] 全页宽 ≤ 170mm — 除 Fig5/6/7 略超需缩放，余均满足 ✓
- [ ] 缩至 170mm 后字体 ≥ 7pt — 全部通过验证 ✓
- [ ] 线条粗细 0.25-1pt — 脚本中 linewidth=0.5~2.0，大部分在范围内
- [ ] Alt Text 已包含 — 已添加至 manuscript_revised_v10.docx 图注中 ✓
- [ ] 图片作为独立文件（非嵌入 docx）— 已提取到此目录 ✓
- [ ] 图注（Figure Legends）在主稿件中 — 在 manuscript_revised_v10.docx 中 ✓
- [ ] 补充文件符合 GigaScience 命名规范（含 "Supplementary Material"）— 已更新 ✓
- [ ] 稿件正文引用使用 "Additional file X" 格式 — 已更新 ✓
- [ ] "List of additional files" 章节已添加 — 已更新 ✓
