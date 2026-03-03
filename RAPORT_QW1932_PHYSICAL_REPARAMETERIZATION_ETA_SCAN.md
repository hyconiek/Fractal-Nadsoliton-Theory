# RAPORT QW-1932: PHYSICAL REPARAMETERIZATION ETA SCAN

- Data UTC: 2026-03-03T05:36:58.534106+00:00
- Verdict: **PHYSICAL_REPARAMETERIZATION_STRICT_PASS**
- strict_pass_count: 6

## Selected
- eta: 2.8
- omega: 0.373414
- phi: -1.310234
- beta: 0.615938
- rel_loss_vs_eta1: -0.9005
- primary corr/gain ratios: 1.0338/2.3426
- stress corr/gain ratios: 1.1587/1.3579
- delta_bic_vs_eta1: -2318.7898

## Strict Flags (selected)
- beta_le_1p20: True
- rel_loss_le_0p35: True
- primary_corr_ratio_ge_0p95: True
- primary_gain_ratio_ge_1p00: True
- stress_corr_ratio_ge_0p95: True
- stress_gain_ratio_ge_1p00: True
- omega_interior: True
- eta_not_grid_edge: True
- delta_bic_le_10: True

## Required Next Step
- INTEGRATE_REPARAM_BRANCH_IN_STAGE_B_GATE_AND_RETEST_IDENTIFIABILITY

## Artifacts
- JSON: `report_qw1932_physical_reparameterization_eta_scan.json`
