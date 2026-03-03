# RAPORT QW-2044: POINTWISE MICRO Z_BETA(d) / DELTA_ETA(d) DERIVATION

- Data UTC: 2026-03-03T21:27:27.667377+00:00
- Verdict: **POINTWISE_MICRO_DERIVATION_FAIL**
- Readiness: **POINTWISE_SUPPORT_NOT_ESTABLISHED**
- pass_count: 4/7

## Target (Refrozen QW-2039)
- beta_target: 0.920000
- eta_target: 1.800000
- Z_beta_target: 92.000000
- delta_eta_target: 0.800000

## Global Micro Estimates
- n_profiles: 54
- n_windows_total: 79
- beta median/CI95: 0.372941 / [0.056120, 41.274830]
- eta median/CI95: 2.937500 / [1.000000, 3.000000]
- Z_beta median: 37.294107
- delta_eta median: 1.937500
- median window rmse: 0.032243

## Pointwise Coverage
- n_bins: 4
- beta target in CI95 fraction: 1.0000
- eta target in CI95 fraction: 0.2500
- joint target in CI95 fraction: 0.2500

## Flags
- enough_pointwise_bins_ge_6: False
- global_beta_target_inside_ci95: True
- global_eta_target_inside_ci95: True
- bin_beta_target_coverage_ge_0p50: True
- bin_eta_target_coverage_ge_0p50: False
- bin_joint_target_coverage_ge_0p35: False
- median_window_rmse_le_0p12: True

## Representative Bins (top-8)
- d=7 n=8 | beta=3.3630 [0.2595,35.0899] | eta=2.3938 [1.0525,2.8413] | rmse=0.0108 | joint=True
- d=6 n=6 | beta=0.8146 [0.0072,18.0110] | eta=2.9312 [2.2359,3.0000] | rmse=0.0024 | joint=False
- d=4 n=18 | beta=3.2410 [0.0976,25.0105] | eta=2.9688 [2.2228,3.0000] | rmse=0.0426 | joint=False
- d=2 n=34 | beta=0.2036 [0.1035,2.1801] | eta=3.0000 [2.4853,3.0000] | rmse=0.0593 | joint=False

## Required Next Step
- IMPROVE_POINTWISE_IDENTIFIABILITY_WITH_PHASE_CONDITIONED_WINDOWS

## Artifacts
- JSON: `report_qw2044_pointwise_micro_zbeta_deltaeta_derivation.json`
