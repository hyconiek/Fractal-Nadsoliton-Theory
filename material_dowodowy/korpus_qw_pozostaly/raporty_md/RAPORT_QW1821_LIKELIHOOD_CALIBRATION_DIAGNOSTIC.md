# RAPORT QW-1821: LIKELIHOOD CALIBRATION DIAGNOSTIC

- Data UTC: 2026-03-02T21:52:50.771096+00:00
- Miscalibration index: 0.554
- Verdict: **LIKELIHOOD_MISSPECIFICATION_SIGNAL_STRONG**
- Recommendation: **Switch primary gate metric to predictive score robust to variance tails (e.g., CRPS / quantile loss / Student-t with calibrated df), then retest.**

## QW-1817
- n_rep: 14
- mean delta LL: 5.4094
- std delta LL: 4.0894
- P(delta LL > 0): 0.929
- mean RMSE gain: 0.258834
- P(RMSE gain > 0): 1.000
- discordance rate: 0.071
- corr(delta LL, RMSE gain): 0.894

## QW-1818
- n_rep: 14
- mean delta LL: 7.7083
- std delta LL: 5.6562
- P(delta LL > 0): 0.857
- mean RMSE gain: 0.316021
- P(RMSE gain > 0): 1.000
- discordance rate: 0.143
- corr(delta LL, RMSE gain): 0.903

## Flags
- ll_unstable: True
- rmse_stable_positive: True
- metric_discordance_present: True

## Artifacts
- JSON: `report_qw1821_likelihood_calibration_diagnostic.json`
