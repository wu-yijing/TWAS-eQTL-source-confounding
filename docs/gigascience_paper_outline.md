# GigaScience 投稿论文大纲

## 题目（建议）

**eQTL source confounding systematically biases TWAS cross-population replication: a quantitative evaluation in HOTAIR binding proteins across 60 genes and three diabetic complications**

或：

**Quantifying the impact of eQTL source selection on TWAS cross-population replication: lessons from HOTAIR binding proteins**

---

## 绪论 (Introduction)

### 段落1: TWAS 作为基因发现工具的成功与局限
- TWAS 已成为将 GWAS 关联定位到基因的重要方法
- 但 TWAS 跨人群复制率普遍偏低（引用已有文献）
- 方法学层面已识别的混杂来源：LD 差异、MAF 差异、eQTL 组织错配
- **本文聚焦的被低估的混杂来源**：eQTL 权重来源的选择（小样本组织特异性 vs 大样本全血）

### 段落2: HOTAIR 案例的引入
- HOTAIR 在糖尿病微血管并发症中的角色（简短，1-2句）
- 前期研究：RNA pull-down 鉴定 30 个 HOTAIR 结合蛋白
- 初始 TWAS 分析（方案1）在 GTEx v8 权重下报告了 15 个 FDR 显著基因-并发症配对，RNH1 为 Top 1（Z=+8.50）
- 但这套初始结果是否稳健？跨人群复制是否成功？

### 段落3: 本研究的核心问题与设计
- 核心科学问题：**eQTL 权重来源的选择在多大程度上影响 TWAS 跨人群复制的结论？**
- 研究设计：以 HOTAIR 为模型系统，系统比较 GTEx v8（组织特异性 MASHR, N≈600）与 eQTLGen（全血, N=31,684）两种权重来源
- 分析框架：104 个基因（30 候选 + 44 非候选 + 30 T2DM 对照，有重叠）× 2 种 eQTL 权重 × 3 个表型（DR/DN/DPN）
  - 非候选组原始 54 个样本特异性蛋白，经 Unused≥2.0 置信度过滤后保留 44 个高置信度基因
- 关键评估指标：方向一致性、Spearman 相关性、富集率变化、协变量匹配后富集稳定性

### 段落4: 主要发现预览
- Spearman ρ = 0.29（所有基因-表型配对）
- 方向一致性仅 59.5%（全基因集）/ 37.5%（候选基因原始19个分析基因子集）
- 富集率从 GTEx 的 40.7%（候选组）下降至 eQTLGen 的 0-7%（候选组 FDR 显著数为 0）
- RNH1 从 Z=+8.50（GTEx）衰减至 Z=+1.18, P=0.237（eQTLGen）
- **结论：eQTL 来源混杂是 TWAS 跨人群复制中的一个主要但被低估的偏倚来源**

---

## 方法 (Methods)

### 2.1 研究设计与分析框架
- 六步分析流程（Figure 1 流程图的文字描述）
- S1 → S2: GTEx baseline TWAS → eQTLGen 敏感性分析
- S3 → S4: Mahalanobis 匹配 → 匹配后富集验证
- S5: Pull-down vs 文献分层分析
- S6: 跨人群验证 + 异质性评估

### 2.2 数据来源
- **FinnGen R13**：DR (15,353/398,626), DN (6,043/396,166), DPN (4,145/393,402)
- **GTEx v8**：Nerve_Tibial (N=532) + Whole_Blood (N=670) MASHR multi-tissue weights
- **eQTLGen**：全血 cis-eQTL (N=31,684)，用于敏感性分析
- **1000 Genomes EUR**：LD 参考面板

### 2.3 S-PrediXcan 分析
- 软件版本、参数（引用 Barbeira et al. 2018）
- GTEx: MASHR 模型，多组织 Stouffer 合并 + ACAT-O 敏感性验证
- eQTLGen: 使用 eQTLGen Z-score 元分析权重
- 多重检验校正：FDR (BH, q<0.05), Bonferroni, 每个分析域分别校正

### 2.4 非候选蛋白置信度过滤（方法学改进）

这是方案2引入的关键改进之一：

- **输入**：54 个非候选蛋白（原始 Covariate Matrix 标注为 `39_NonCandidate_HOTAIR`，但实际含 54 个基因——该标签为历史命名遗留问题）
- **标准**：Unused 评分 ≥ 2.0（蛋白鉴定置信度阈值，替代无区分力的 Peptides95≥1 标准）
- **输出**：44 个高置信度非候选蛋白，排除 10 个低置信度蛋白（RPS10P5, RPL7A, SERPINH1, RPL17, H2AJ, RPL13, TPM4, TUBB4B, VAT1, H3-7）
- 该过滤确保了 Mahalanobis 匹配的对照基因质量，排除了可能因鉴定不可靠而引入噪声的低置信度蛋白

### 2.5 两组内部对照 + 第三组 T2DM 对照
- 候选组：30 HOTAIR RNA pull-down 蛋白（Unused≥2.0 + CRAPome 过滤）
- 非候选组：44 高置信度非候选蛋白（原始 54 经 Unused≥2.0 过滤后保留，排除 10 个低置信度蛋白）
- T2DM 对照组：30 已发表 T2DM GWAS 关联基因
- 分组详细依据（Supplementary Table S1）

### 2.6 方向一致性与 Spearman 相关性分析
- 逐基因-表型配对比较 GTEx Z 值与 eQTLGen Z 值的符号一致性
- 整体 Spearman 秩相关 ρ（合并三表型）
- 分组分表型报告方向一致率 + 二项符号检验

### 2.7 富集分析
- FDR 显著率 = (组内 FDR 显著基因数) / (组内可测试基因数)
- 富集比较：Fisher 精确检验（候选 vs 非候选 / 候选 vs T2DM）
- 分三层：GTEx baseline → eQTLGen → Mahalanobis 匹配后

### 2.8 Mahalanobis 协变量匹配（方法学核心）
- 目标：为 30 个候选基因从全基因组基因池中匹配协变量相近的对照
- 匹配维度：log10(基因长度) + GC% + eQTL SNP 数
- 方法：matchit R 包，nearest neighbor, Mahalanobis 距离, ratio=1
- 平衡诊断：SMD 比较（Love plot, Supplementary Figure S2）
- 匹配后：在匹配配对集上重复富集分析

### 2.9 Pull-down 来源分层分析
- 将 30 候选基因分为两组：RNA pull-down 实验来源 vs 文献/数据库来源
- 分别计算两组的 eQTLGen 方向一致性和 Z 值分布

### 2.10 跨人群复制与异质性评估
- UK Biobank 两独立 DR GWAS 复制
- 补充外部验证：Cai et al. 2026 UKB T2D-DR, Sakaue et al. 2021 DN
- Cochran's Q, I², 95% 预测区间

### 2.11 统计软件与可用性
- 所有分析均使用 Python 3.13 + R 4.5.2
- 关键包：S-PrediXcan (v1.0), matchit R, scipy, statsmodels, scikit-learn
- 完整代码已公开在 GitHub: [URL]

---

## 结果 (Results)

### 3.1 eQTL 来源选择对 TWAS Z 值的系统影响
**【核心发现展示】**

- **方向一致性仅约 60%**：
  - 全基因集（合并三表型）：59.5% 方向一致
  - 候选基因（19 个 GTEx 可分析基因）：仅 37.5%
  
- **Spearman 相关性极弱**：
  - 整体 ρ = 0.289 (P<0.001)
  - 分表型：DR ρ=0.31, DN ρ=0.27, DPN ρ=0.28（Figure 2）

- **前导信号 RNH1 完全衰减**：
  - GTEx: Z=+8.50 → eQTLGen: Z=+1.18, P=0.237
  - 所有核心基因 Z 值变化（Table 1 或 Figure 3）

- **基因方向翻转普遍**：CKAP4, HSP90AB1, HSP90B1, H2AFV, XRCC5

### 3.2 富集分析结论因 eQTL 来源选择而反转
**【方法学核心发现】**

- **GTEx 框架下**：候选组 FDR 显著率 40.7%，高于非候选组 33.3% 和 T2DM 对照 21.1%（Fisher 未显著）
- **eQTLGen 框架下**：候选组 FDR 显著率 0-7%（三表型合并 0%），非候选组 38-76%，T2DM 50-79%
- **关键反转**：eQTLGen 权重下非候选组富集率**高于**候选组——初始 TWAS 富集结论完全不可复制
- **Fisher 精确检验**：因候选组无 FDR 显著基因而无法计算

### 3.3 Mahalanobis 匹配验证
- 30 候选基因全部成功匹配（ratio=1）
- 匹配后协变量 SMD 均 <0.1（Figure S2, Love Plot）
- 匹配后的富集分析：候选组仍无 FDR 显著基因
- **结论**：富集差异的消失并非由基因长度/GC/eQTL SNP 数的协变量差异驱动，而是 eQTL 来源本身的系统效应

### 3.4 GTEx 与 eQTLGen 的富集率全面对比（Table 2）

| 表型 | 权重来源 | 候选组 FDR率 | 非候选组 FDR率 | T2DM组 FDR率 |
|------|---------|-------------|---------------|-------------|
| DR | GTEx | 40.7% | 33.3% | 21.1% |
| DR | eQTLGen | 0% | 38.1-50% | 50-78.6% |
| DN | GTEx | ... | ... | ... |
| DN | eQTLGen | ... | ... | ... |
| DPN | GTEx | ... | ... | ... |
| DPN | eQTLGen | 0% | 42.9% | 50% |

→ **可视化**：Figure 5 富集率对比条形图

### 3.5 Pull-down vs 文献分层分析
- 30 候选 → 15 pull-down 来源 + 15 文献/数据库来源
- 两组之间无系统性方向一致性差异（Figure S3）
- **含义**：候选来源不影响 eQTL 混杂的方向——混杂效应是全局性的

### 3.6 跨人群复制：独立验证结果
- UK Biobank DR：8/8 未复制（方向一致 87.5%, P=0.035；但全部 P>0.05）
- RNH1 补充验证：方向一致但效应量衰减 3.7-15 倍
- DR 三队列 Meta 分析：FE Z=+3.37 (P=7.4×10⁻⁴), I²=95.0%, 95% PI [-4.00, +10.68]
- 与 eQTLGen 敏感性分析一致：跨人群复制失败的部分原因是 UKB 使用 eQTLGen 权重而 FinnGen 使用 GTEx 权重

### 3.7 综合证据总结（Table 3）
| 假设 | GTEx 证据 | eQTLGen 证据 | Mahalanobis 证据 | 复制证据 | 综合判定 |
|------|-----------|-------------|-----------------|---------|---------|
| RNH1 跨并发症候选 | 强 (Z=+8.50) | 无效 (Z=+1.18) | 无效 | 方向一致但不可复制 | 需要重新评估 |
| 候选组富集 | 40.7% vs 33.3% | 0% vs 38-76% | 0% | — | 不成立 |
| 跨人群复制 | 87.5% 方向一致 | — | — | I²=95% | 不可复制 |

---

## 讨论 (Discussion)

### 4.1 TWAS 中 eQTL 来源混杂的系统性
- 本研究的核心方法学发现：在分析框架中 eQTL 来源的选择（GTEx vs eQTLGen）系统性改变 TWAS 的定性结论
- 方向一致性 ~60% 意味着仅凭一种 eQTL 权重的 TWAS 发现约 40% 的信号是 eQTL 来源依赖性的
- **建议**：未来 TWAS 研究必须报告至少两种独立 eQTL 权重来源的敏感性分析；仅使用单一 eQTL 来源的 TWAS 发现应标注"eQTL 来源依赖性"

### 4.2 为什么 eQTLGen 与 GTEx 差异如此之大？
- **样本量差异**：eQTLGen (N=31,684) vs GTEx (N≈600)，统计功效差异巨大 → 理论上 eQTLGen 更可靠
- **组织特异性**：GTEx 提供神经+血液，eQTLGen 仅提供全血，但大样本量弥补了部分组织错配
- **权重估计方法**：MASHR 贝叶斯降维 vs 单变量线性元分析 → 在 LD 模式和效应量异质性下的表现差异
- **真正的混杂来源于**：TWAS 的 instrumental variable（eQTL SNP）选择因 eQTL 来源不同而不同，导致同一基因在不同权重系统中被赋予不同的遗传代理变量

### 4.3 对 HOTAIR 候选基因生物学意义的重新评估
- 简短讨论（2-3 句）：原始候选筛选策略的价值（在 GTEx 中确实看到候选 > T2DM 的趋势）和局限（eQTL 混杂无法通过生物学先验校正）
- 对 RNH1 的定位：从"Top 1 跨并发症候选"修正为"需 eQTL 来源敏感性分析验证的假设性靶点"
- Pull-down vs 文献来源的无差异进一步支持了 eQTL 混杂是全局性问题

### 4.4 与已有 TWAS 方法学文献的比较
- 已有文献关注：组织错配、SNP 稀疏选择偏倚、共定位功效、Stouffer 独立性假设
- 本文新增层面：**eQTL 权重来源选择本身的系统性效应**——不同于组织错配（关注不同组织）或方法差异（关注统计模型），本研究关注的是"同一组织但不同数据集/样本量的 eQTL 权重"对结论的影响

### 4.5 本研究的边界条件
- FinnGen 芬兰隔离人群的特异性——始祖效应的放大效应
- 仅涵盖糖尿病微血管并发症三种表型
- 仅测试两种 eQTL 权重来源（GTEx v8 + eQTLGen）
- 缺乏独立第三组织（如视网膜）作为金标准对照
- 未纳入 Biobank Japan, All of Us 等非欧洲人群

### 4.6 对未来 TWAS 研究的建议（方法学贡献的落地）
1. **强制性敏感性分析**：所有 TWAS 研究应报告至少两种独立 eQTL 来源的结果
2. **方向一致性报告**：报告跨 eQTL 来源的 Spearman ρ 和方向一致率
3. **分层报告**：按 eQTL 来源分层报告富集分析
4. **代码和数据可用性**（呼应 GigaScience 定位）：所有分析脚本和数据标准化输出

### 4.7 结论
- HOTAIR 案例证明：eQTL 来源选择是 TWAS 跨人群复制中一个被低估的关键混杂因素
- 方向一致性 ~60% 意味着仅基于单一 eQTL 权重的 TWAS 发现中约 40% 可能是 eQTL 来源依赖性的
- 建议将双源 eQTL 敏感性分析纳入 TWAS 研究的报告标准

---

## 图表清单 (Figures & Tables)

### Main Figures (6)
| 编号 | 标题 | 数据来源 | 建议类型 |
|------|------|---------|---------|
| Figure 1 | 研究设计与分析流程图 | 设计文档 | 流程图/示意图 |
| Figure 2 | GTEx vs eQTLGen TWAS Z 值散点图（三表型并排） | eqtlgen_vs_gtex_comparison.csv | 散点图 |
| Figure 3 | 核心基因双源 Z 值对比瀑布图 | eqtlgen_spredixcan_results.csv | 条形图/瀑布图 |
| Figure 4 | 三组富集率对比：GTEx vs eQTLGen | enrichment_comparison.csv | 分组条形图 |
| Figure 5 | Mahalanobis 匹配前后协变量平衡（Love Plot） | mahalanobis_matched_pairs.csv + covariate matrix | Love Plot |
| Figure 6 | RNH1 跨人群/跨 eQTL 来源效应衰减综合图 | 多来源合并 | 综合图/森林图+衰减图 |

### Main Tables (3-4)
| 编号 | 标题 |
|------|------|
| Table 1 | GTEx vs eQTLGen 核心基因 Z 值与方向一致性对比（涵盖 RNH1, RPL10A, SERPINH1 等核心基因） |
| Table 2 | 三组 × 三表型 × 两种 eQTL 来源的 FDR 富集率完整矩阵 |
| Table 3 | 综合证据评估表（最终结论导向） |

### Supplementary Materials
| 编号 | 内容 |
|------|------|
| Figure S1 | Stouffer vs ACAT-O 敏感性验证散点图 |
| Figure S2 | Pull-down vs 文献分层方向一致性对比 |
| Figure S3 | Mahalanobis 匹配 27/30 配对完整对比 |
| Figure S4 | DR/DN/DPN 三表型 |Z| 密度分布（eQTLGen） |
| Table S1 | 完整协变量矩阵 |
| Table S2 | 原始 GTEx TWAS 完整结果表（方案1基线） |
| Table S3 | eQTLGen TWAS 完整结果表 |
| Table S4 | Mahalanobis 匹配配对完整表 |
| Table S5 | 富集分析完整结果 |
| Table S6 | 跨人群复制完整结果 + 异质性统计 |

---

## 关键数据流（从输入到输出的映射）

```
数据输入层：
  ├─ FinnGen R13 GWAS (DR/DN/DPN Z-score)
  ├─ GTEx v8 MASHR weights (Nerve_Tibial + Whole_Blood)
  ├─ eQTLGen whole blood cis-eQTL (Z-score meta-analysis)
  ├─ 1000G EUR LD reference
  ├─ covariate matrix (gene length, GC%, source annotation)
  │
S-PrediXcan层：
  ├─ GTEx-based TWAS (baseline)
  └─ eQTLGen-based TWAS (sensitivity)
       │
分析层：
  ├─ Direction consistency analysis
  ├─ Spearman correlation
  ├─ Enrichment analysis (per group × per phenotype × per eQTL source)
  ├─ Mahalanobis matching (matchit R)
  ├─ Post-matching enrichment validation
  └─ Pull-down vs literature stratification
       │
输出层（定量指标）：
  ├─ Spearman ρ = 0.29 (all pairs)
  ├─ Direction consistency = 59.5% (all) / 37.5% (candidate)
  ├─ FDR rate: GTEx 40.7% → eQTLGen 0% (candidate)
  ├─ RNH1: Z=+8.50 → Z=+1.18 (P=0.237)
  ├─ Mahalanobis SMD < 0.1 (balanced)
  └─ Cross-population I² = 95.0%
       │
方法学结论：
  └─ eQTL source confounding is a systematic bias in TWAS
     └→ Future TWAS must report dual-source sensitivity analysis
```

---

## 建议的Cover Letter关键信息

- **为何选 GigaScience**：本文核心贡献是方法学发现 + 开放式数据/代码 + FAIR 原则对齐
- **不选 JTM（原方案1目标期刊）的原因**：原稿受限于单一 eQTL 来源、单一人群，GigaScience 更适合重新定位后的方法学 + 数据驱动型论文
- **作者需确认**：
  - 所有数据来自公共资源，符合使用条款
  - 数据和代码已在 GigaDB + GitHub 公开
  - 所有作者无利益冲突

---

*大纲生成日期: 2026-06-27 | 目标期刊: GigaScience | 文章类型: Research Article*
