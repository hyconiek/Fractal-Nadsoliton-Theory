# RAPORT QW-1818: ROBUST SEQUENCE OOS

- Data UTC: 2026-03-02T21:48:10.296034+00:00
- OOS cohort size: 90
- Mean windows per pair: 15.59
- Features: f_mean, f_std, f_slope, f_quad, f_spread, f_autoc1, f_switch
- Replications: 14
- Mean test delta LL (M2E-M2): 7.7083
- P(test LL M2E>M2): 0.857
- Mean RMSE gain (M2-M2E): 0.316021
- P(RMSE gain>0): 1.000
- Std delta LL (QW-1817 -> QW-1818): 4.089 -> 5.656
- Std reduction ratio: -0.383
- Verdict: **ROBUST_SEQUENCE_OOS_PARTIAL**

## Pass Flags
- test_ll_gain: True
- test_rmse_gain: True
- dispersion_reduction_vs_qw1817: False

## Artifacts
- JSON: `report_qw1818_robust_sequence_oos.json`
