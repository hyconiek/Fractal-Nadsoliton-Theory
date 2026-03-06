# RAPORT QW-2018: V2 TRIAD PERMUTATION CONVERGENCE AUDIT

- Data UTC: 2026-03-03T19:15:47.884501+00:00
- Verdict: **V2_TRIAD_PERMUTATION_CONVERGENCE_ROBUST_PASS**
- Required next step: `PROCEED_TO_TRIAD_STAGE_WITH_V2_PACKAGE`

## Holdout Base Metrics
- pearson_corr: 0.053107
- spearman_corr: 0.030768
- rmse_gain_ratio: 0.000275

## Permutation Grid
- n_perm=1000: p_corr_med=0.01049 [q10=0.00799, q90=0.01389], p_gain_med=0.00699 [q10=0.00310, q90=0.01179], frac_corr<=0.01=0.500, frac_gain<=0.01=0.833
- n_perm=3000: p_corr_med=0.00833 [q10=0.00633, q90=0.01193], p_gain_med=0.00616 [q10=0.00367, q90=0.00796], frac_corr<=0.01=0.583, frac_gain<=0.01=1.000
- n_perm=5000: p_corr_med=0.00910 [q10=0.00736, q90=0.01054], p_gain_med=0.00590 [q10=0.00502, q90=0.00678], frac_corr<=0.01=0.833, frac_gain<=0.01=1.000
- n_perm=10000: p_corr_med=0.00925 [q10=0.00813, q90=0.01047], p_gain_med=0.00565 [q10=0.00485, q90=0.00648], frac_corr<=0.01=0.750, frac_gain<=0.01=1.000
- n_perm=20000: p_corr_med=0.00937 [q10=0.00830, q90=0.00979], p_gain_med=0.00590 [q10=0.00536, q90=0.00644], frac_corr<=0.01=0.917, frac_gain<=0.01=1.000

## Robust Flags
- robust_corr_p_le_0p01: True
- robust_gain_p_le_0p01: True

## Artifacts
- JSON: `report_qw2018_v2_triad_permutation_convergence_audit.json`
