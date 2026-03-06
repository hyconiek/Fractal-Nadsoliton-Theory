# RAPORT QW-2045: PHASE-CONDITIONED POINTWISE MICRO DERIVATION

- Data UTC: 2026-03-03T21:29:42.874057+00:00
- Verdict: **PHASE_CONDITIONED_POINTWISE_DERIVATION_PARTIAL**
- Readiness: **PARTIAL_IDENTIFIABILITY_REPAIR**
- pass_count: 6/8

## Target (Refrozen QW-2039)
- beta_target: 0.920000
- eta_target: 1.800000

## Global Estimates
- n_profiles: 66
- n_rows_total: 81
- beta median/CI95: 0.333205 / [0.056034, 37.175785]
- eta median/CI95: 2.675000 / [1.450000, 3.000000]
- rmse median: 0.033877
- phase_min median: 0.593487

## Pointwise Coverage
- n_bins: 2
- beta target in CI95 fraction: 1.0000
- eta target in CI95 fraction: 1.0000
- joint target in CI95 fraction: 1.0000

## Flags
- enough_pointwise_bins_ge_6: False
- global_beta_target_inside_ci95: True
- global_eta_target_inside_ci95: True
- bin_beta_target_coverage_ge_0p50: True
- bin_eta_target_coverage_ge_0p50: True
- bin_joint_target_coverage_ge_0p35: True
- median_window_rmse_le_0p10: True
- phase_condition_strength_ge_0p75: False

## Representative Bins (top-8)
- d=2 n=56 | beta=0.3414 [0.1398,40.7554] | eta=2.6250 [1.4484,3.0000] | rmse=0.0268 | phase=0.595 | joint=True
- d=4 n=25 | beta=0.2843 [0.0437,30.8342] | eta=3.0000 [1.4800,3.0000] | rmse=0.0928 | phase=0.554 | joint=True

## Required Next Step
- RUN_MICRO_STAGEC_INTERSECTION_GATE

## Artifacts
- JSON: `report_qw2045_phase_conditioned_pointwise_micro_derivation.json`
