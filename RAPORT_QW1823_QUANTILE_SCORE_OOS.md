# RAPORT QW-1823: QUANTILE SCORE OOS

- Data UTC: 2026-03-02T21:57:48.039718+00:00
- OOS cohort size: 90
- Mean windows per pair: 15.59
- Features: f_mean, f_std, f_slope, f_quad, f_spread, f_autoc1, f_switch
- Taus: [0.1, 0.5, 0.9]
- Replications: 14
- Mean quantile gain (M2-M2E): 0.062399
- P(quantile gain>0): 1.000
- Mean MAE gain (M2-M2E): 0.350754
- P(MAE gain>0): 1.000
- Std quantile gain: 0.025727
- Verdict: **QUANTILE_SCORE_OOS_SUPPORTED**

## Pass Flags
- quantile_gain: True
- mae_gain: True
- dispersion_control: True

## Artifacts
- JSON: `report_qw1823_quantile_score_oos.json`
