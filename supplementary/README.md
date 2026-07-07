# GigaScience 补充材料 — 文件目录

## 目录结构

```
上传论文supplyment/
├── README.md                        ← 本文件
├── figures/                         ← 补充图片（PDF + PNG + 图注，更新后命名）
│   ├── Figure_S1_Supplementary_Material.pdf    ← 30对 Mahalanobis 匹配配对图
│   ├── Figure_S2_Supplementary_Material.pdf    ← Pull-down vs Literature 分层富集图
│   ├── Figure_S3_Supplementary_Material.pdf    ← Stouffer vs ACAT-O 多组织整合敏感性
│   ├── Figure_S4_Supplementary_Material.pdf    ← |Z| 密度分布（eQTLGen）
│   ├── Figure_S5_Supplementary_Material.pdf    ← eQTL SNP 数分布
│   ├── Figure_S6_Supplementary_Material.pdf    ← 方向一致性敏感性分析
│   ├── Figure_S7_Supplementary_Material.pdf    ← 分组 |Z| 密度分布
│   ├── Figure_S8_Supplementary_Material.pdf    ← 所有表型 FDR 富集汇总
│   └── 对应 .png 和 .txt（图注）文件
├── tables/                          ← 补充表格（CSV + XLSX + 表注，更新后命名）
│   ├── Table_S1_Supplementary_Material.csv    ← 104 基因分组列表
│   ├── Table_S2_Supplementary_Material.csv    ← GTEx baseline TWAS 结果
│   ├── Table_S3_Supplementary_Material.csv    ← 基因协变量矩阵
│   ├── Table_S4_Supplementary_Material.csv    ← Mahalanobis 匹配配对
│   ├── Table_S5_Supplementary_Material.csv    ← 跨人群复制结果
│   ├── Table_1_gtex_vs_eqtlgen_comparison.csv  ← 主文 Table 1
│   ├── Table_2_fdr_enrichment_rates.csv       ← 主文 Table 2
│   ├── Table_3_integrated_evidence_assessment.csv ← 主文 Table 3
│   └── 对应 .xlsx 和 .txt（表注）文件
└── data/                            ← GigaDB 数据文件
    ├── gigadb_metadata_form.csv     ← GigaDB 元数据表单
    ├── README.md                    ← GigaDB 数据说明
    ├── eqtlgen_spredixcan_results.csv
    ├── eqtlgen_vs_gtex_comparison.csv
    ├── covariate_matrix.csv
    ├── mahalanobis_matched_pairs.csv
    ├── enrichment_comparison.csv
    ├── candidate_comparison_DR.csv
    ├── layer_analysis.csv
    └── viz_z_distribution.csv
```

---

## 与稿件的对应关系（GigaScience Additional file 格式）

### Additional Files 编号对应

| Additional File # | 稿件中引用 | 实际文件 | 说明 |
|:-----------------:|-----------|---------|------|
| Additional file 1 | Figure S1 | `figures/Figure_S1_Supplementary_Material.pdf` | 30对匹配配对图 |
| Additional file 2 | Figure S2 | `figures/Figure_S2_Supplementary_Material.pdf` | Pull-down vs Literature 分层 |
| Additional file 3 | Figure S3 | `figures/Figure_S3_Supplementary_Material.pdf` | ACAT-O 敏感性 |
| Additional file 4 | Figure S4 | `figures/Figure_S4_Supplementary_Material.pdf` | |Z| 密度分布 |
| Additional file 5 | Figure S5 | `figures/Figure_S5_Supplementary_Material.pdf` | eQTL SNP 数分布 |
| Additional file 6 | Figure S6 | `figures/Figure_S6_Supplementary_Material.pdf` | 方向一致性敏感性 |
| Additional file 7 | Figure S7 | `figures/Figure_S7_Supplementary_Material.pdf` | 分组 |Z| 密度分布 |
| Additional file 8 | Figure S8 | `figures/Figure_S8_Supplementary_Material.pdf` | 所有表型 FDR 富集汇总 |
| Additional file 9 | Table S1 | `tables/Table_S1_Supplementary_Material.csv` | 104 基因分组列表 |
| Additional file 10 | Table S2 | `tables/Table_S2_Supplementary_Material.csv` | GTEx baseline TWAS |
| Additional file 11 | Table S3 | `tables/Table_S3_Supplementary_Material.csv` | 基因协变量矩阵 |
| Additional file 12 | Table S4 | `tables/Table_S4_Supplementary_Material.csv` | Mahalanobis 匹配配对 |
| Additional file 13 | Table S5 | `tables/Table_S5_Supplementary_Material.csv` | 跨人群复制结果 |

---

## 文件大小汇总

| 类别 | 数量 | 总大小 |
|------|:----:|:------:|
| 补充图片 PDF | 8 | 250 KB |
| 补充图片 PNG | 8 | 1.6 MB |
| 补充表格 CSV | 8 (含主文Table) | 55 KB |
| 图例/表注 TXT | 13 | 5 KB |
| GigaDB 数据 | 10 | 52 KB |
| **合计** | **47** | **~2 MB** |

> 所有文件均 ≤20 MB，符合 GigaScience 附加文件大小要求。

---

**备注：**
- `Table_1`/`Table_2`/`Table_3` 为主文中的表格，以独立 CSV 文件存放备用
- 稿件中所有引用已更新为 `[see Additional file X: Figure/Table SX]` 格式
- GigaDB DOI 待插入 Data Availability 部分
- 文件命名符合 GigaScience "supplementary material" 命名规范

*Updated on 2026-07-07 for GigaScience submission*
