# -*- coding: utf-8 -*-
"""Fig. 4 — axis-resolved source-discordance decomposition, 49-gene common universe
(147 gene-phenotype pairs; the same universe as the tissue-only arm of Table 3b).

PROVENANCE OF EVERY PLOTTED VALUE  (read before reusing this script)
--------------------------------------------------------------------------
(a) Reproduced from source data and re-verified 2026-09-14:
      panel-only  rho = 0.355   direction consistency = 67.3%
      tissue-only rho = 0.422   direction consistency = 67.3%
      dual-arm direction consistency = 59.2%  ; d(direction consistency) = -8.2 pp
(b) Transcribed verbatim from the archived analysis record (2026-09-11 validation
    snapshot, consistent with the manuscript) because they could NOT be recomputed
    exactly from the archived inputs:
      dual-arm rho = 0.287 with its 95% CI, and the delta-rho / delta-P values of panel (b).

Why (b) is not recomputed here: ACAT-O combines two-sided tissue-level P-values and
therefore carries no intrinsic direction. The sign convention used in the original run was
not recorded in the archived code; it is now documented in Methods 2.4 (the sign is taken
from the GTEx Whole_Blood Z-score of the same gene-phenotype pair). Applying that
convention reproduces the direction-consistency figures exactly, but yields dual-arm
rho = 0.364 rather than 0.287, and no combination of archived P-value sources tested with
that sign convention reproduces 0.287. The manuscript reports the archived value and this
figure plots it verbatim, so figure and text remain consistent; re-deriving the dual arm
from raw inputs is an open item (see 出图脚本硬编码审计报告_20260914.md).

Before publishing a new version of this package, either (i) re-derive the dual arm and
replace the literals, or (ii) keep this provenance block so the transcription is explicit.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8

OUT = r'E:\workbuddy\BMC Genomics投稿资料\定稿图集_Fig1-8_20260914'
os.makedirs(OUT, exist_ok=True)

# ===== plotted values: see the PROVENANCE block in the module docstring =====
# 49-gene common universe, n = 147 pairs
arm_names = ['Panel-only\nresource / sample-size\n67.3% consistent',
             'Tissue-only\ntissue context\n67.3% consistent',
             'Dual\nboth axes\n59.2% consistent']
# (a) recomputed: panel-only and tissue-only rho ; (b) transcribed: dual-arm rho
rho = np.array([0.355, 0.422, 0.287])
ci_lo = np.array([0.204, 0.279, 0.131])
ci_hi = np.array([0.488, 0.547, 0.429])
# (a) recomputed
cons = np.array([67.3, 67.3, 59.2])
ci_err = np.vstack([rho - ci_lo, ci_hi - rho])

# (b) Forest (panel b): transcribed from the archived analysis record
labels_b = ['dual − panel-only', 'dual − tissue-only']
delta_rho = np.array([-0.067, -0.135])
rho_lo = np.array([-0.196, -0.375])
rho_hi = np.array([+0.039, +0.081])
p_rho = np.array([0.252, 0.238])
delta_dc = np.array([-8.2, -8.2])
p_dc = np.array([0.018, 0.201])

C_PANEL, C_TISSUE, C_DUAL = '#1f77b4', '#d62728', '#9467bd'
TXT = '#222222'

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.4, 3.35), dpi=300,
                                 gridspec_kw={'width_ratios': [1.0, 1.30], 'wspace': 0.42})

# ---------------- (a) 柱状图 ----------------
x = np.arange(3)
ax_a.bar(x, rho, width=0.62, color=[C_PANEL, C_TISSUE, C_DUAL],
         edgecolor='black', linewidth=0.6, zorder=3)
ax_a.errorbar(x, rho, yerr=ci_err, fmt='none', ecolor='black',
              capsize=4, capthick=0.9, elinewidth=0.9, zorder=4)

for xi, r, hi, lo in zip(x, rho, ci_hi, ci_lo):
    ax_a.text(xi, hi + 0.026, r'$\rho$ = %.3f' % r, ha='center', va='bottom',
              fontsize=8.2, fontweight='bold', color=TXT, zorder=5)
    ax_a.text(xi, hi + 0.088, '95%% CI %.3f–%.3f' % (lo, hi), ha='center', va='bottom',
              fontsize=6.0, color='#555555', zorder=5)
for xi, c in zip(x, cons):
    pass

ax_a.set_xticks(x)
ax_a.set_xticklabels(arm_names, fontsize=6.8)
ax_a.set_xlim(-0.62, 2.62)
ax_a.set_ylim(0, 0.80)
ax_a.set_yticks(np.arange(0, 0.61, 0.1))
ax_a.set_ylabel('Spearman $\\rho$ between sources', fontsize=8.5)
ax_a.tick_params(axis='y', labelsize=7.5)
ax_a.tick_params(axis='x', length=0)
ax_a.text(0.02, 0.97, '(a)', transform=ax_a.transAxes, ha='left', va='top',
          fontsize=10, fontweight='bold', color='black', zorder=10)
ax_a.text(0.98, 0.985, '49-gene common universe\nn = 147 pairs', transform=ax_a.transAxes,
          ha='right', va='top', fontsize=6.4, color='#444444', style='italic', zorder=10)
for sp in ('top', 'right'):
    ax_a.spines[sp].set_visible(False)
for sp in ('left', 'bottom'):
    ax_a.spines[sp].set_color('black'); ax_a.spines[sp].set_linewidth(0.8)

# ---------------- (b) Forest ----------------
y = np.array([1, 0])
ax_b.axvline(0, color=C_TISSUE, linestyle=(0, (5, 3)), linewidth=1.0, zorder=2)
for yi, lo, hi in zip(y, rho_lo, rho_hi):
    ax_b.plot([lo, hi], [yi, yi], color='black', linewidth=1.0, zorder=3, solid_capstyle='butt')
    ax_b.plot([lo, lo], [yi - 0.10, yi + 0.10], color='black', linewidth=1.0, zorder=3)
    ax_b.plot([hi, hi], [yi - 0.10, yi + 0.10], color='black', linewidth=1.0, zorder=3)
ax_b.scatter(delta_rho, y, s=40, color='black', edgecolor='black', linewidth=0.5, zorder=4)

for yi, dr, pr, dd, pd_, lo, hi in zip(y, delta_rho, p_rho, delta_dc, p_dc, rho_lo, rho_hi):
    ax_b.text(0.145, yi + 0.46,
              r'$\Delta\rho$ = %+.3f  [%+.3f, %+.3f], P = %.3f' % (dr, lo, hi, pr),
              ha='right', va='center', fontsize=6.6, color=TXT,
              bbox=dict(boxstyle='round,pad=0.30', facecolor='white', edgecolor='none', alpha=0.95),
              zorder=10)
    ax_b.text(0.145, yi + 0.20,
              r'$\Delta$ consistency = %+.1f pp, P = %.3f' % (dd, pd_),
              ha='right', va='center', fontsize=6.6, color='#555555',
              bbox=dict(boxstyle='round,pad=0.30', facecolor='white', edgecolor='none', alpha=0.95),
              zorder=10)

ax_b.set_yticks(y)
ax_b.set_yticklabels(labels_b, fontsize=7.6)
ax_b.set_ylim(-0.65, 1.75)
ax_b.set_xlim(-0.55, 0.15)
ax_b.set_xticks(np.arange(-0.5, 0.11, 0.1))
ax_b.set_xlabel('Difference in $\\rho$ (dual minus single-axis)', fontsize=8.5)
ax_b.tick_params(axis='x', labelsize=7.5)
ax_b.tick_params(axis='y', length=0)
ax_b.text(0.02, 0.97, '(b)', transform=ax_b.transAxes, ha='left', va='top',
          fontsize=10, fontweight='bold', color='black', zorder=10)
for sp in ('top', 'right'):
    ax_b.spines[sp].set_visible(False)
for sp in ('left', 'bottom'):
    ax_b.spines[sp].set_color('black'); ax_b.spines[sp].set_linewidth(0.8)

fig.subplots_adjust(left=0.085, right=0.985, bottom=0.235, top=0.93)

# 自检：文字不得溢出
_r = fig.canvas.get_renderer(); W, H = fig.canvas.get_width_height(); bad = []
for _ax in fig.axes:
    for _t in _ax.texts:
        bb = _t.get_window_extent(renderer=_r)
        if bb.x0 < 3 or bb.x1 > W - 3 or bb.y0 < 3 or bb.y1 > H - 3:
            bad.append('"%s" x=%.0f~%.0f y=%.0f~%.0f' % (_t.get_text()[:34], bb.x0, bb.x1, bb.y0, bb.y1))
print('[自检] ' + ('全部文字在画布内' if not bad else '溢出警告:\n  ' + '\n  '.join(bad)))

for ext, kw in [('png', dict(dpi=600, facecolor='white')), ('pdf', dict(facecolor='white'))]:
    p = os.path.join(OUT, 'Fig4.%s' % ext)
    plt.savefig(p, **kw)
    print('已生成:', p)
plt.close()

print()
print('数值核对: (a) ρ = 0.355 / 0.422 / 0.287 ; 方向一致 = 67.3 / 67.3 / 59.2')
print('          (b) Δρ = -0.067 / -0.135 ; Δ方向一致 = -8.2 / -8.2 pp')
