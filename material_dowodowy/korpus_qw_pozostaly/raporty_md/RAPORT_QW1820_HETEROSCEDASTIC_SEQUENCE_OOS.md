# RAPORT QW-1820: HETEROSCEDASTIC SEQUENCE OOS

- Data UTC: 2026-03-02T21:51:40.983972+00:00
- OOS cohort size: 90
- Mean windows per pair: 15.59
- Features: f_mean, f_std, f_slope, f_quad, f_spread, f_autoc1, f_switch
- Replications: 14
- Mean delta LL (hetero-homo): -1.8629
- P(hetero>homo): 0.357
- Mean delta LL (hetero-M2): 5.4567
- P(hetero>M2): 0.786
- Std delta LL vs M2 (QW-1817 -> hetero): 4.089 -> 5.780
- Std reduction ratio: -0.413
- Verdict: **HETEROSCEDASTIC_SEQUENCE_OOS_WEAK**

## Pass Flags
- hetero_gain_vs_homo: False
- hetero_gain_vs_m2: False
- dispersion_reduction_vs_qw1817: False

## Artifacts
- JSON: `report_qw1820_heteroscedastic_sequence_oos.json`
