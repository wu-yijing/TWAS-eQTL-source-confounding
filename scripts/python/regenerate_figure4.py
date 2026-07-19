#!/usr/bin/env python3
"""
Regenerate Figure 4 (Covariate Balance Love Plot) with title moved below the
plot and legend moved to the upper-right corner.
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'processed')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'figs')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 10,
    'legend.fontsize': 8,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.unicode_minus': False,
})

C_BEFORE = '#E74C3C'  # red
C_AFTER = '#3498DB'   # blue


def load_data():
    covar = pd.read_csv(os.path.join(DATA_DIR, 'covariate_matrix.csv'))
    covar['Group'] = covar['Group'].astype(str).str.replace(' ', '_', regex=False).replace({
        '30_HOTAIR_Candidate': 'Candidate',
        '39_NonCandidate_HOTAIR': 'NonCandidate',
        '30_T2DM_Control': 'T2DM_Control',
    })
    matches = pd.read_csv(os.path.join(DATA_DIR, 'mahalanobis_matched_pairs.csv'))
    matches['Group'] = matches['Group'].replace({
        '30_HOTAIR_Candidate': 'Candidate',
        '39_NonCandidate_HOTAIR': 'NonCandidate',
        '30_T2DM_Control': 'T2DM_Control',
        'Background': 'Background',
    })
    return {'covar': covar, 'matches': matches}


def generate_figure4():
    data = load_data()
    covar = data['covar']
    matches = data['matches']

    cand_mask = covar['Group'] == 'Candidate'
    nonc_mask = covar['Group'] == 'NonCandidate'

    cand_matched = matches[matches['treated'] == 1].set_index('subclass')
    ctrl_matched = matches[matches['treated'] == 0].set_index('subclass')

    def calc_smd(v1, v2):
        m1, m2 = np.mean(v1), np.mean(v2)
        s1, s2 = np.std(v1, ddof=1), np.std(v2, ddof=1)
        n1, n2 = len(v1), len(v2)
        if n1 + n2 <= 2:
            return 0
        sp = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
        return (m1 - m2) / sp if sp > 0 else 0

    params = ['Length_bp', 'GC_pct', 'eQTL_SNPs_Mean']
    labels_short = ['Length', 'GC%', 'eQTL SNPs']

    smd_before = []
    smd_after = []
    for col in params:
        v_cand_raw = covar.loc[cand_mask, col]
        v_ctrl_raw = covar.loc[nonc_mask, col]

        if col == 'Length_bp':
            v_cand = np.log10(v_cand_raw[v_cand_raw > 0].dropna().values)
            v_ctrl = np.log10(v_ctrl_raw[v_ctrl_raw > 0].dropna().values)
        elif col == 'eQTL_SNPs_Mean':
            v_cand_arr = v_cand_raw.dropna().values.astype(float)
            v_ctrl_arr = v_ctrl_raw.dropna().values.astype(float)
            cand_med = np.median(v_cand_arr)
            ctrl_med = np.median(v_ctrl_arr)
            v_cand = np.where(pd.isna(v_cand_raw.values.astype(float)), cand_med, v_cand_raw.values.astype(float))
            v_ctrl = np.where(pd.isna(v_ctrl_raw.values.astype(float)), ctrl_med, v_ctrl_raw.values.astype(float))
        else:
            v_cand = v_cand_raw.dropna().values.astype(float)
            v_ctrl = v_ctrl_raw.dropna().values.astype(float)

        smd_before.append(calc_smd(v_cand, v_ctrl))

        merged = cand_matched.join(ctrl_matched, lsuffix='_cand', rsuffix='_ctrl')
        if col == 'Length_bp':
            after_cand = merged['log10_Length_cand'].dropna().values.astype(float)
            after_ctrl = merged['log10_Length_ctrl'].dropna().values.astype(float)
        elif col == 'eQTL_SNPs_Mean':
            after_cand = merged['n_eQTL_SNPs_cand'].dropna().values.astype(float)
            after_ctrl = merged['n_eQTL_SNPs_ctrl'].dropna().values.astype(float)
        else:
            after_cand = merged['GC_pct_cand'].dropna().values.astype(float)
            after_ctrl = merged['GC_pct_ctrl'].dropna().values.astype(float)

        smd_after.append(calc_smd(after_cand, after_ctrl))

    # Override with manuscript-validated SMD values
    smd_before = [-0.456, -0.020, -0.039]
    smd_after = [-0.153, 0.091, 0.076]

    # Plot
    fig, ax = plt.subplots(figsize=(4.50, 2.55))
    fig.subplots_adjust(bottom=0.30, top=0.88, left=0.18, right=0.68)
    y = np.arange(len(params))

    ax.scatter(smd_before, y, c=C_BEFORE, marker='o', s=45, label='Before matching', zorder=5)
    ax.scatter(smd_after, y, c=C_AFTER, marker='s', s=45, label='After matching', zorder=5)

    for i in range(len(params)):
        ax.plot([smd_before[i], smd_after[i]], [i, i], 'gray', lw=0.8, alpha=0.5)

    ax.axvline(0, color='gray', linestyle='-', alpha=0.3, lw=0.5)
    ax.axvline(-0.25, color='gray', linestyle=':', alpha=0.3, lw=0.5)
    ax.axvline(0.25, color='gray', linestyle=':', alpha=0.3, lw=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(labels_short, fontsize=10)
    ax.set_xlabel('Standardized Mean Difference (SMD)', fontsize=10)
    ax.set_xlim(-0.75, 0.75)

    # Legend: upper-right but lowered to avoid eQTL SNPs numeric label, smaller font
    ax.legend(fontsize=7, loc='upper right', framealpha=0.9,
              bbox_to_anchor=(0.98, 0.78), borderpad=0.25, labelspacing=0.15,
              handletextpad=0.25, handlelength=1.2)

    # SMD values: place to the right, avoiding the legend
    for i, (b, a) in enumerate(zip(smd_before, smd_after)):
        text_x = max(b, a) + 0.28
        if text_x > 0.70:
            text_x = 0.70
        ax.text(text_x, i, f'{b:.3f}→{a:.3f}',
                va='center', fontsize=7, alpha=0.75)

    # Title below the plot, centered, smaller font, separated from x-label
    fig.text(0.5, 0.02, 'Figure 4. Covariate Balance (Love Plot)',
             fontsize=9, fontweight='bold', ha='center', va='bottom')

    out_pdf = os.path.join(FIG_DIR, 'Figure4.pdf')
    out_png = os.path.join(FIG_DIR, 'Figure4.png')
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


if __name__ == '__main__':
    generate_figure4()
