# RAPORT QW-1867: BETA-AUGMENTED OBSERVABLES BENCHMARK

- Data UTC: 2026-03-03T00:06:06.776900+00:00
- Verdict: **BETA_AUGMENTATION_SUPPORTED**
- n_rep: 170

## Arms
- A_BASELINE: paired=False | beta RMSE=0.0180 | beta tol=0.676 | corr_ob=-0.075
- B_AUGMENTED: paired=False | beta RMSE=0.0109 | beta tol=0.847 | corr_ob=0.042
- C_AUGMENTED_PAIRED: paired=True | beta RMSE=0.0152 | beta tol=0.765 | corr_ob=-0.311

## Ranking (vs baseline)
- B_AUGMENTED: score=1.416, beta_rmse_factor=1.654, beta_tol_gain=0.171, |corr| reduction=0.034
- C_AUGMENTED_PAIRED: score=1.094, beta_rmse_factor=1.187, beta_tol_gain=0.088, |corr| reduction=-0.236
- A_BASELINE: score=1.000, beta_rmse_factor=1.000, beta_tol_gain=0.000, |corr| reduction=0.000

## Best Arm
- B_AUGMENTED

## Artifacts
- JSON: `report_qw1867_beta_augmented_observables_benchmark.json`
