# RAPORT QW-1817: SEQUENCE OOS VALIDATION

- Data UTC: 2026-03-02T21:45:55.653877+00:00
- OOS cohort size: 90
- Mean windows per pair: 15.59
- Features: f_mean, f_std, f_slope, f_quad, f_spread, f_autoc1, f_switch
- Replications: 14
- Mean test delta LL (M2E-M2): 5.4094
- P(test LL M2E>M2): 0.929
- Mean RMSE gain (M2-M2E): 0.258834
- P(RMSE gain>0): 1.000
- Std test delta LL: 4.089
- Verdict: **SEQUENCE_OOS_VALIDATION_PARTIAL**

## Pass Flags
- test_ll_gain: True
- test_rmse_gain: True
- dispersion_control: False

## Artifacts
- JSON: `report_qw1817_sequence_oos_validation.json`
