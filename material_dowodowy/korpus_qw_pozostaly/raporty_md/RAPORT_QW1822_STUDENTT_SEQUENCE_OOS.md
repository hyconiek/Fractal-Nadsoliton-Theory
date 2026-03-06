# RAPORT QW-1822: STUDENT-T SEQUENCE OOS

- Data UTC: 2026-03-02T21:55:32.354247+00:00
- OOS cohort size: 90
- Mean windows per pair: 15.59
- Features: f_mean, f_std, f_slope, f_quad, f_spread, f_autoc1, f_switch
- Replications: 14
- Mean delta t-LL (M2E-M2): 7.1490
- P(delta t-LL > 0): 0.786
- P(M2E_t > flat_t): 0.786
- Mean nu(M2E): 25.714
- Std delta LL (QW-1817 -> Student-t): 4.089 -> 6.114
- Std reduction ratio: -0.495
- Verdict: **STUDENTT_SEQUENCE_OOS_WEAK**

## Pass Flags
- studentt_gain_vs_m2: False
- studentt_gain_vs_flat: False
- dispersion_reduction_vs_qw1817: False

## Artifacts
- JSON: `report_qw1822_studentt_sequence_oos.json`
