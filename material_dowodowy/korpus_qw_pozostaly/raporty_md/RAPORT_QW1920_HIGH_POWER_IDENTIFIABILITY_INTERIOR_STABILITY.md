# RAPORT QW-1920: HIGH-POWER IDENTIFIABILITY INTERIOR STABILITY

- Data UTC: 2026-03-03T04:25:16.402994+00:00
- Verdict: **HIGH_POWER_IDENTIFIABILITY_DATA_LIMITED_INTERIOR_NOT_ESTABLISHED**
- Required next step: `DESIGN_ORTHOGONAL_BETA_OBSERVABLE_ON_BLIND_EXTERNAL_WITH_INTERVENTION_PROTOCOL`

## Real Arm (extended profiles)
- n_profiles: 42
- n_points_per_profile: 24
- optimum: omega=0.056250, phi=0.942500, beta=2.000000
- unique_modes: 3
- hessian_cond: 1.986e+16
- boundary beta_high fraction: 1.000

## Synthetic Power Arm
- n_rep: 180
- sigma_noise: 0.0500
- joint_hit_rate: 1.000
- boundary_estimate_rate: 0.033
- median_abs_error omega/phi/beta: 0.0024 / 0.0284 / 0.0025

## Flags
- estimator_power_strong: True
- interior_stable_real: False

## Artifacts
- JSON: `report_qw1920_high_power_identifiability_interior_stability.json`
