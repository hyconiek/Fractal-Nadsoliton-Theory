# RAPORT QW-2048: SPECTRAL PHASE-LOCKED POINTWISE DERIVATION

- Data UTC: 2026-03-03T23:37:02.367015+00:00
- Verdict: **SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION_PASS**
- Readiness: **POINTWISE_IDENTIFIABILITY_REPAIRED**
- pass_count: 8/8

## Target (Refrozen QW-2039)
- beta_target: 0.920000
- eta_target: 1.800000

## Global Estimates
- n_rows_total: 342
- beta median/CI95: 1.147396 / [0.000054, 106.682605]
- eta median/CI95: 2.125000 / [1.000000, 3.000000]
- rmse median: 0.004257
- phase_min median: 0.852165
- omega_phase median (q10,q90): 0.277107 (0.176524, 0.414192)
- phase_r2 median: 0.088173
- n_informative median: 10.00

## Pointwise Coverage
- n_bins: 17
- beta target in CI95 fraction: 0.8235
- eta target in CI95 fraction: 0.9412
- joint target in CI95 fraction: 0.8235

## Flags
- enough_pointwise_bins_ge_6: True
- global_beta_target_inside_ci95: True
- global_eta_target_inside_ci95: True
- bin_beta_target_coverage_ge_0p50: True
- bin_eta_target_coverage_ge_0p50: True
- bin_joint_target_coverage_ge_0p35: True
- median_window_rmse_le_0p10: True
- phase_condition_strength_ge_0p75: True

## Representative Bins (top-10)
- d=17 n=11 | beta=0.1133 [0.0000,2.9079] | eta=2.6125 [1.0000,3.0000] | rmse=0.0001 | phase=0.848 | joint=True
- d=18 n=13 | beta=0.1408 [0.0003,2.3093] | eta=2.1875 [1.0000,3.0000] | rmse=0.0001 | phase=0.880 | joint=True
- d=14 n=12 | beta=0.2851 [0.0023,113.1427] | eta=1.1562 [1.0000,2.8522] | rmse=0.0002 | phase=0.847 | joint=True
- d=15 n=10 | beta=0.4152 [0.0000,5.5777] | eta=1.8938 [1.0000,3.0000] | rmse=0.0003 | phase=0.864 | joint=True
- d=16 n=10 | beta=2.1994 [0.0246,46.8576] | eta=2.1812 [1.2128,2.9747] | rmse=0.0003 | phase=0.889 | joint=True
- d=2 n=85 | beta=0.4354 [0.1585,63.2250] | eta=2.2625 [1.2400,3.0000] | rmse=0.0005 | phase=0.889 | joint=True
- d=13 n=13 | beta=1.0976 [0.0052,31.2584] | eta=1.4875 [1.0000,2.9663] | rmse=0.0009 | phase=0.881 | joint=True
- d=3 n=29 | beta=0.2724 [0.0601,184.1459] | eta=2.5375 [1.4563,3.0000] | rmse=0.0038 | phase=0.845 | joint=True
- d=12 n=27 | beta=2.4765 [0.0008,72.5430] | eta=1.8000 [1.0000,3.0000] | rmse=0.0427 | phase=0.842 | joint=True
- d=11 n=16 | beta=6.3624 [0.0024,762.6814] | eta=1.8375 [1.0703,3.0000] | rmse=0.0674 | phase=0.857 | joint=True

## Required Next Step
- RUN_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE

## Artifacts
- JSON: `report_qw2048_spectral_phase_locked_pointwise_derivation.json`
