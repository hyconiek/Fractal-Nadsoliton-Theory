# RAPORT QW-1941: TRIPLE-SECTOR SHARED COMPLEXITY (AIC/BIC)

- Data UTC: 2026-03-03T06:19:30.790610+00:00
- Verdict: **TRIPLE_SECTOR_NOT_CLOSED_UNDER_SHARED_COMPLEXITY_CONTROL**
- Chosen model: **M1_ONE_SHARED_LAMBDA**

## M0 Strict Derived
- total_loss: 11.9842
- all_pass: False
- AIC/BIC proxy: 139.0812/139.0812

## M1 One Shared Lambda (best)
- lambda: 0.7000
- total_loss: 8.4045
- all_pass: False
- AIC/BIC proxy: 121.2110/122.5432

## Comparison
- delta_aic (M1-M0): -17.8702
- delta_bic (M1-M0): -16.5380

## Required Next Step
- REDESIGN_SHARED_FLAVOR_MASS_LINK_UNDER_SINGLE_KERNEL_CONSTRAINT

## Artifacts
- JSON: `report_qw1941_triple_sector_shared_complexity_aic_bic.json`
