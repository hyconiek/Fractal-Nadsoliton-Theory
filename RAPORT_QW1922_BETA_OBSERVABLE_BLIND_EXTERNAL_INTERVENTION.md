# RAPORT QW-1922: BETA OBSERVABLE BLIND EXTERNAL INTERVENTION

- Data UTC: 2026-03-03T04:28:13.691278+00:00
- Selected observable: `B7_local_resid_std`
- Verdict: **BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PASS**

## Primary (holdout)
- effect_beta: 0.9907
- effect_omega: 0.0435
- contrast: 0.9472
- contrast_boot_q05: 0.8338
- all_pass: True

## Stress (holdout)
- effect_beta: 2.4068
- effect_omega: 0.1888
- contrast: 2.2180
- contrast_boot_q05: 2.0730
- all_pass: True

## Global Flags
- primary_all_pass: True
- stress_all_pass: True

## Required Next Step
- INTEGRATE_BETA_OBSERVABLE_AS_EXPLICIT_CONSTRAINT_IN_TRIAD_DERIVATION

## Artifacts
- JSON: `report_qw1922_beta_observable_blind_external_intervention.json`
