#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M6(b)/M6(d) sensitivity analyses: N_e-weighted re-merge of the RNH1 DR cross-cohort meta-analysis.

Background (Methods 1.10 of the manuscript): the primary merge combines cohort Z-scores as an
UNWEIGHTED arithmetic mean (each cohort enters with unit variance, SE = 1). This script repeats
the merge with sqrt(N_e) weights as a pre-specified sensitivity analysis, and additionally
reports the descriptive three-study merge that re-includes the underpowered UK Biobank dataset
GCST90043640 that was excluded from the primary estimate.

Framework
---------
Under a common standardized effect, an S-PrediXcan Z-score scales approximately with
sqrt(N_e):  Z_i ~ N(theta * sqrt(N_e,i), 1).
Hence the per-sqrt(N_e)-unit effect is  b_i = Z_i / sqrt(N_e,i)  with precision N_e,i, and
fixed-effect inverse-variance pooling of the b_i is algebraically identical to a weighted
Stouffer combination of the Z-scores with weights sqrt(N_e,i):

    Z_meta = sum_i Z_i * sqrt(N_e,i) / sqrt(sum_i N_e,i)

Cochran's Q is evaluated on the b scale:  Q = sum_i N_e,i * (b_i - b_hat)^2.
Random effects use the DerSimonian-Laird moment estimator for tau^2 on the b scale.

Inputs (all previously published in the manuscript / archived outputs)
----------------------------------------------------------------------
  FinnGen R13 DR   : Z = +13.318 (eQTLGen weights; data/processed/eqtlgen_spredixcan_harmonized_results.csv), N_e = 59,134
  UKB ieu-b-4803   : Z = +0.95  (eQTLGen weights; Xue et al. 2022),          N_e = 54,209 = 4*14147*322390/336537
  UKB GCST90043640 : Z = +0.57  (GTEx v8 MASHR Nerve_Tibial weights),        N_e =  1,231 = 4*308*456040/456348
                     (excluded from the primary estimate; underpowered, different weight source)

Outputs
-------
  stdout / data/processed/m6_ne_weighted_sensitivity_results.txt
"""
import math

STUDIES = [
    ("FinnGen R13 DR (eQTLGen weights)", 13.318, 59134),
    ("UKB ieu-b-4803 (eQTLGen weights)", 0.95, 54209),
    ("UKB GCST90043640 (GTEx v8 MASHR Nerve_Tibial weights)", 0.57, 1231),
]


def pnorm_two_sided(z):
    return math.erfc(abs(z) / math.sqrt(2))


def weighted_stouffer(zs, nes):
    num = sum(z * math.sqrt(n) for z, n in zip(zs, nes))
    return num / math.sqrt(sum(nes))


def ivw_b(zs, nes):
    """Fixed-effect IVW of per-sqrt(N_e)-unit effects; returns (b_hat, Q)."""
    bs = [z / math.sqrt(n) for z, n in zip(zs, nes)]
    b_hat = sum(n * b for n, b in zip(nes, bs)) / sum(nes)
    Q = sum(n * (b - b_hat) ** 2 for n, b in zip(nes, bs))
    return b_hat, Q


def dl_random_effects(zs, nes):
    """DerSimonian-Laird random effects on the b scale."""
    b_hat, Q = ivw_b(zs, nes)
    k = len(zs)
    tau2 = max(0.0, (Q - (k - 1)) /
               (sum(nes) - sum(n * n for n in nes) / sum(nes)))
    w = [1.0 / (1.0 / n + tau2) for n in nes]
    b_re = sum(wi * (z / math.sqrt(n)) for wi, z, n in zip(w, zs, nes)) / sum(w)
    se_re = math.sqrt(1.0 / sum(w))
    z_re = b_re / se_re
    return b_re, se_re, tau2, z_re, pnorm_two_sided(z_re)


def report(label, zs, nes):
    k = len(zs)
    df = k - 1
    zw = weighted_stouffer(zs, nes)
    b_hat, Q = ivw_b(zs, nes)
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    b_re, se_re, tau2, z_re, p_re = dl_random_effects(zs, nes)
    mean_unw = sum(zs) / k
    lines = [
        "=" * 78, label,
        "  studies: " + "; ".join("Z=%+.3f, N_e=%d" % (z, n) for z, n in zip(zs, nes)),
        "  N_e shares: " + ", ".join("%.1f%%" % (100.0 * n / sum(nes)) for n in nes),
        "  [primary reference] unweighted mean Z = %+.2f" % mean_unw,
        "  [sqrt(N_e)-weighted] pooled Z = %+.2f  (two-sided P = %.2g)" % (zw, pnorm_two_sided(zw)),
        "  [sqrt(N_e)-weighted] b-scale FE pooled b = %.5f" % b_hat,
        "  [heterogeneity] Cochran Q = %.1f (df=%d), I^2 = %.1f%%" % (Q, df, I2),
        "  [random effects] DL tau^2 = %.3g (b scale); pooled b = %.4f (SE %.4f); "
        "test Z = %.2f, P = %.3f" % (tau2, b_re, se_re, z_re, p_re),
    ]
    print("\n".join(lines))
    return lines


def main():
    lines = []
    # M6(b): primary k=2 set (FinnGen + ieu-b-4803), sqrt(N_e)-weighted re-merge
    lines += report("M6(b)  k=2 primary set, sqrt(N_e)-weighted re-merge",
                    [13.318, 0.95], [59134, 54209])
    # M6(d): descriptive k=3 merge including GCST90043640
    lines += report("M6(d)  k=3 all-studies merge incl. GCST90043640",
                    [13.318, 0.95, 0.57], [59134, 54209, 1231])
    # unweighted k=3 reference
    zs = [13.318, 0.95, 0.57]
    m3 = sum(zs) / 3
    Q3 = sum((z - m3) ** 2 for z in zs)
    lines += ["=" * 78,
              "REFERENCE  k=3 unweighted mean: Z = %+.2f ; Q = %.1f (df=2) ; I^2 = %.1f%%"
              % (m3, Q3, (Q3 - 2) / Q3 * 100)]
    text = "\n".join(lines) + "\n"
    print()
    import os
    out = os.path.join(os.path.dirname(__file__), "..", "..",
                       "data", "processed", "m6_ne_weighted_sensitivity_results.txt")
    with open(out, "w", encoding="utf-8") as f:
        f.write(text)
    print("Saved:", os.path.normpath(out))


if __name__ == "__main__":
    main()
