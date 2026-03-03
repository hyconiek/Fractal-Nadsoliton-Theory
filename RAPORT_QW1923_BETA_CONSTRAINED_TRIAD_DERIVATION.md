# RAPORT QW-1923: BETA-CONSTRAINED TRIAD DERIVATION

- Data UTC: 2026-03-03T04:30:03.182655+00:00
- Verdict: **BETA_CONSTRAINED_TRIAD_DERIVATION_TRADEOFF_HIGH**
- Required next step: `TUNE_CONSTRAINT_WEIGHT_AND_RETEST_TRANSFER`

## Beta Prior (from QW-1922)
- effect_obs_primary: 0.9907
- effect_obs_stress: 2.4068
- beta_target: 0.2837
- beta_interval: [0.1180, 0.4870]
- beta_scale: 0.1845

## Fit Comparison
- unconstrained beta: 2.0000
- constrained beta: 0.6050
- unconstrained objective: 104.0602
- constrained objective: 253.9028
- relative objective increase: 1.4400
- beta_boundary_relief: True
- fit_preserved: False

## Artifacts
- JSON: `report_qw1923_beta_constrained_triad_derivation.json`
