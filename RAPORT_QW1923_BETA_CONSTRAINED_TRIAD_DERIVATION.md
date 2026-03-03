# RAPORT QW-1923: BETA-CONSTRAINED TRIAD DERIVATION

- Data UTC: 2026-03-03T19:06:33.181379+00:00
- Verdict: **BETA_CONSTRAINED_TRIAD_DERIVATION_TRADEOFF_HIGH**
- Required next step: `TUNE_CONSTRAINT_WEIGHT_AND_RETEST_TRANSFER`

## Beta Prior (from QW-1922)
- effect_obs_primary: 1.4071
- effect_obs_stress: 2.4068
- beta_target: 0.3515
- beta_interval: [0.1240, 0.4810]
- beta_scale: 0.1785

## Fit Comparison
- unconstrained beta: 2.0000
- constrained beta: 0.6325
- unconstrained objective: 104.0602
- constrained objective: 241.1974
- relative objective increase: 1.3179
- beta_boundary_relief: True
- fit_preserved: False

## Artifacts
- JSON: `report_qw1923_beta_constrained_triad_derivation.json`
