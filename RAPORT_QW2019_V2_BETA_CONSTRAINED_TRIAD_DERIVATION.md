# RAPORT QW-2019: BETA-CONSTRAINED TRIAD DERIVATION

- Data UTC: 2026-03-03T19:16:36.547461+00:00
- Verdict: **BETA_CONSTRAINED_TRIAD_DERIVATION_TRADEOFF_HIGH**
- Required next step: `TUNE_CONSTRAINT_WEIGHT_AND_RETEST_TRANSFER`

## Beta Prior (from QW-1922)
- effect_obs_primary: 1.4381
- effect_obs_stress: 2.4068
- beta_target: 0.1632
- beta_interval: [0.0832, 0.2432]
- beta_scale: 0.0800

## Fit Comparison
- unconstrained beta: 2.0000
- constrained beta: 0.3141
- unconstrained objective: 104.0602
- constrained objective: 363.8723
- relative objective increase: 2.4967
- beta_boundary_relief: True
- fit_preserved: False

## Artifacts
- JSON: `report_qw2019_v2_beta_constrained_triad_derivation.json`
