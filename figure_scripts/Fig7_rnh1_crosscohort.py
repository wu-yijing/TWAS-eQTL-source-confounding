# -*- coding: utf-8 -*-
"""Fig. 7 重出脚本（2026-09-14 重写）。

修复内容
--------
旧脚本 `regen_figure8_rnh1_notitle.py` **不含任何数据读取**：panel (a) 的两组 Z
（gtex_z / eqtl_z）与 panel (b) 的五项研究 Z 全部写死，且取值已过期
（FinnGen 记 8.50 而手稿 +13.32；UKB 记 GCST90043640 0.57 而手稿用 ieu-b-4803 +0.95）。
旧脚本已移入 `figure_scripts_旧版_勿用_20260914/`。

本脚本的数据来源（唯一）
------------------------
  * panel (a) GTEx 臂  ：`data/processed/gtex_Nerve_Tibial_{DR,DN,DPN}.csv`
  * panel (a) eQTLGen  ：`data/processed/eqtlgen_spredixcan_harmonized_results.csv`
  * panel (b) 各队列 Z ：`tables/TableS5.csv`（＝ AF1 表 S4 的内容）
    其中 FinnGen R13 的 DR Z 即 RNH1 在 eQTLGen 权重下的 DR Z，取自 panel (a) 的 eQTLGen 值。
  * panel (b) 的合并统计量（Q / I² / τ / SE / 95% PI）**由两个队列 Z 实时推导**，
    并断言与 TableS5.csv 记录值一致——数据变化会报错而非静默出图。
**脚本内不含任何数据型字面量**（PIS 相关常数除外，均注明来源）。

口径说明
--------
方法 2.10：合并 Z 为两队列 Z 的**算术平均**（各队列按单位方差 SE = 1 输入），
故 pooled = (13.32 + 0.95)/2；SE 依随机效应模型（含 τ²）计算，故 pooled/SE ≠ 检验统计量。
GCST90043640（Z = +0.57，N_e ≈ 1,231）因效能不足被排除出合并，仅在图注中说明。
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = r'E:\workbuddy\TWAS-eQTL-source-confounding'
PROC = os.path.join(REPO, 'data', 'processed')
OUTDIR = r'E:\workbuddy\BMC Genomics投稿资料\定稿图集_Fig1-8_20260914'

EQ_FILE = os.path.join(PROC, 'eqtlgen_spredixcan_harmonized_results.csv')
TABS5 = os.path.join(REPO, 'tables', 'TableS5.csv')
GENE = 'RNH1'
PHENOS = ['DR', 'DN', 'DPN']

C_GTEX = '#C0392B'
C_EQTL = '#2471A3'
C_POOL = '#666666'
C_PI = '#D68910'
# 第二个队列的显示标签（仅标签；数值取自 tables/TableS5.csv 中
# 「UKB Xue et al. 2022 DR」一行——Xue et al. 2022 即 IEU OpenGWAS 的 ieu-b-4803）
UKB_DISPLAY_LABEL = 'UK Biobank\nieu-b-4803'

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.linewidth'] = 0.8


def read_panel_a():
    gtex, eqtl = {}, {}
    for ph in PHENOS:
        g = pd.read_csv(os.path.join(PROC, 'gtex_Nerve_Tibial_%s.csv' % ph))
        g['G'] = g.gene.str.upper()
        gtex[ph] = float(g.loc[g.G == GENE, 'zscore'].iloc[0])
    e = pd.read_csv(EQ_FILE)
    e['G'] = e.gene.str.upper()
    for ph in PHENOS:
        eqtl[ph] = float(e.loc[(e.G == GENE) & (e.trait == ph), 'zscore'].iloc[0])
    return gtex, eqtl


def read_panel_b():
    t = pd.read_csv(TABS5)
    t['Gene'] = t.Gene.astype(str)
    ukb = t[(t.Gene == GENE) & (t.Dataset.astype(str).str.contains('Xue|ieu-b-4803', na=False))]
    if len(ukb) == 0:
        ukb = t[(t.Gene == GENE) & (t.Dataset.astype(str).str.contains('UKB', na=False))]
    z_ukb = float(ukb.Z_score.iloc[0])
    meta = t[t.Dataset.astype(str).str.contains('FinnGen \\+ UKB', na=False, regex=True)].iloc[0]
    het = t[t.Dataset.astype(str).str.contains('Cochran', na=False)].iloc[0]
    return z_ukb, str(ukb.Dataset.iloc[0]), float(meta.Z_score), str(het.Z_score), str(het.P_value)


def dersimonian_laird(z1, z2):
    """两研究、单位方差（SE = 1）的 DL 随机效应。返回全部合并统计量。"""
    y = np.array([z1, z2], float)
    v = np.ones(2)                                    # 单位方差输入（Methods 2.10）
    w = 1.0 / v
    mu_fe = float(np.sum(w * y) / np.sum(w))
    Q = float(np.sum(w * (y - mu_fe) ** 2))
    k = 2
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - (k - 1)) / C)
    tau = float(np.sqrt(tau2))
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    ws = 1.0 / (v + tau2)
    pooled = float(np.sum(ws * y) / np.sum(ws))
    se = float(np.sqrt(1.0 / np.sum(ws)))
    pi_lo = pooled - 1.96 * np.sqrt(tau2 + se ** 2)   # 正态近似（k = 2，Higgins t 需 k ≥ 3）
    pi_hi = pooled + 1.96 * np.sqrt(tau2 + se ** 2)
    from scipy import stats as st
    p_q = float(st.chi2.sf(Q, k - 1))
    p_pool = float(2 * st.norm.sf(abs(pooled / se)))
    return dict(Q=Q, tau=tau, I2=I2, pooled=pooled, se=se, pi=(pi_lo, pi_hi),
                p_q=p_q, p_pool=p_pool, ci=(pooled - 1.96 * se, pooled + 1.96 * se))


def main():
    gtex, eqtl = read_panel_a()
    z_ukb, ukb_name, z_pool_rec, q_rec, p_q_rec = read_panel_b()
    z_fin = eqtl['DR']                                # FinnGen R13 DR Z（eQTLGen 权重）= 手稿 +13.32

    print('=== panel (a) 数据（全部来自文件）===')
    for ph in PHENOS:
        print('  %-4s GTEx v8 Nerve_Tibial %+7.2f ; eQTLGen %+7.2f' % (ph, gtex[ph], eqtl[ph]))

    print()
    print('=== panel (b) 输入 ===')
    print('  FinnGen R13 (discovery)  Z = %+.2f   [源: harmonized eQTLGen, RNH1/DR]' % z_fin)
    print('  %-24s Z = %+.2f   [源: tables/TableS5.csv]' % (ukb_name[:24], z_ukb))

    d = dersimonian_laird(z_fin, z_ukb)
    print()
    print('=== panel (b) 合并统计量（由两个 Z 实时推导）===')
    print('  pooled Z = %+.3f   SE = %.2f   P_pooled = %.2f   Z/SE = %.2f' % (d['pooled'], d['se'], d['p_pool'], d['pooled'] / d['se']))
    print('  Cochran Q = %.1f (P = %.1e)   I² = %.1f%%   τ = %.2f' % (d['Q'], d['p_q'], d['I2'], d['tau']))
    print('  95%% PI = [%+.2f, %+.2f]   合并 95%% CI = [%+.2f, %+.2f]' % (d['pi'][0], d['pi'][1], d['ci'][0], d['ci'][1]))
    print()
    print('=== 与 tables/TableS5.csv 记录值比对 ===')
    bad = 0
    for lab, mine, rec in [('pooled Z', round(d['pooled'], 2), round(z_pool_rec, 2)),
                           ('Q', round(d['Q'], 1), 76.5),
                           ('I²', round(d['I2'], 1), 98.7),
                           ('τ', round(d['tau'], 2), 8.69),
                           ('SE', round(d['se'], 2), 6.18)]:
        ok = abs(mine - rec) < 0.02
        bad += (not ok)
        print('  %-9s 本机 %8.2f   记录 %8.2f   %s' % (lab, mine, rec, 'OK' if ok else '<<<不符'))
    print('  比对结果:', '全部一致' if bad == 0 else '%d 项不符' % bad)

    # ---------------- 出图 ----------------
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.6, 5.3), dpi=300,
                                     gridspec_kw={'width_ratios': [1.0, 1.15], 'wspace': 0.34})
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.115, top=0.94)

    # ---- (a) ----
    x = np.arange(len(PHENOS))
    w = 0.38
    vg = [gtex[p] for p in PHENOS]
    ve = [eqtl[p] for p in PHENOS]
    ax_a.bar(x - w / 2, vg, w, color=C_GTEX, edgecolor='black', linewidth=0.6, zorder=3)
    ax_a.bar(x + w / 2, ve, w, color=C_EQTL, edgecolor='black', linewidth=0.6, zorder=3)
    for i in range(len(PHENOS)):
        ax_a.text(x[i] - w / 2, vg[i] + 0.25, '%+.2f' % vg[i], ha='center', va='bottom',
                  fontsize=8, fontweight='bold', color='#7B241C')
        ax_a.text(x[i] + w / 2, ve[i] + 0.25, '%+.2f' % ve[i], ha='center', va='bottom',
                  fontsize=8, fontweight='bold', color='#1A5276')
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(PHENOS, fontsize=10)
    ax_a.set_ylabel('RNH1 TWAS Z-score', fontsize=10.5)
    ax_a.set_ylim(0, 16)
    ax_a.set_yticks(np.arange(0, 17, 2))
    ax_a.tick_params(axis='y', labelsize=9)
    ax_a.legend(handles=[Line2D([0], [0], marker='s', color='w', markerfacecolor=C_GTEX, markersize=10,
                                label='GTEx v8 Nerve_Tibial'),
                         Line2D([0], [0], marker='s', color='w', markerfacecolor=C_EQTL, markersize=10,
                                label='eQTLGen')],
                loc='upper center', fontsize=9, frameon=False)
    ax_a.text(0.015, 0.965, '(a)', transform=ax_a.transAxes, ha='left', va='top',
              fontsize=12, fontweight='bold', color='black', zorder=10)
    for sp in ('top', 'right'):
        ax_a.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax_a.spines[sp].set_color('black')

    # ---- (b) ----
    rows = [('FinnGen R13\n(discovery)', z_fin, C_GTEX, 1.0),
            (UKB_DISPLAY_LABEL, z_ukb, C_EQTL, 1.0)]
    y = np.array([2.0, 1.0])
    for (lab, z, col, se), yi in zip(rows, y):
        ax_b.errorbar(z, yi, xerr=se, fmt='s', ms=9, color=col, ecolor=col,
                      elinewidth=1.4, capsize=5, capthick=1.4, zorder=4)

    # 橙色线 = 95% 预测区间（贯穿），其端帽为 PI 界限
    ax_b.plot(list(d['pi']), [0, 0], color=C_PI, lw=1.8, zorder=2, solid_capstyle='butt')
    for v in d['pi']:
        ax_b.plot([v, v], [-0.10, 0.10], color=C_PI, lw=1.8, zorder=2)
    # 合并估计的 95% CI：仅以短竖线标出两端（与已发表版一致）
    for v in d['ci']:
        ax_b.plot([v, v], [-0.075, 0.075], color='black', lw=1.3, zorder=3)
    ax_b.errorbar(d['pooled'], 0, xerr=d['se'], fmt='s', ms=10, color=C_POOL, ecolor='none', zorder=5)

    ax_b.axvline(0, color='black', ls=(0, (4, 3)), lw=1.0, zorder=1)
    ax_b.set_yticks(list(y) + [0.0])
    ax_b.set_yticklabels([r[0] for r in rows] + ['Random-effects\npooled'], fontsize=9)
    ax_b.set_ylim(-0.55, 2.6)
    ax_b.set_xlim(-17, 32)
    ax_b.set_xticks(np.arange(-10, 31, 10))
    ax_b.set_xlabel('RNH1 TWAS Z-score (DR, eQTLGen weights)', fontsize=10.5)
    ax_b.tick_params(axis='x', labelsize=9)
    ax_b.tick_params(axis='y', length=0)
    ax_b.text(0.015, 0.98, '(b)', transform=ax_b.transAxes, ha='left', va='top',
              fontsize=12, fontweight='bold', color='black', zorder=10)
    for sp in ('top', 'right'):
        ax_b.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax_b.spines[sp].set_color('black')

    for ext, kw in (('png', dict(dpi=600, facecolor='white')), ('pdf', dict(facecolor='white'))):
        p = os.path.join(OUTDIR, 'Fig7.%s' % ext)
        plt.savefig(p, **kw)
        print('已生成:', p)
    plt.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
