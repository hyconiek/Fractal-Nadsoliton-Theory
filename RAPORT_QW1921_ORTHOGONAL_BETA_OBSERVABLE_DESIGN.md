# RAPORT QW-1921: ORTHOGONAL BETA OBSERVABLE DESIGN

- Data UTC: 2026-03-03T04:26:50.668177+00:00
- Verdict: **ORTHOGONAL_BETA_OBSERVABLE_DESIGN_PASS**
- Selected candidate (discovery): `B7_local_resid_std`

## Split
- n_discovery: 985
- n_holdout: 1031

## Selected Candidate Metrics
- discovery corr_abs_with_omega_proxy: 0.0231
- discovery beta_sensitivity: 1.0755
- discovery omega_leakage: 0.0418
- discovery orthogonal_score: 1.0106
- holdout corr_abs_with_omega_proxy: 0.0939
- holdout beta_sensitivity: 0.9907
- holdout omega_leakage: 0.0435
- holdout orthogonal_score: 0.8533

## Pass Flags
- holdout_corr_abs_le_0p20: True
- holdout_beta_sensitivity_ge_0p35: True
- holdout_omega_leakage_le_0p25: True
- holdout_orthogonal_score_positive: True

## Required Next Step
- IMPLEMENT_SELECTED_BETA_OBSERVABLE_IN_BLIND_EXTERNAL_INTERVENTION_RUN

## Artifacts
- JSON: `report_qw1921_orthogonal_beta_observable_design.json`
