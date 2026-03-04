# RAPORT QW-2015: TRUE EXTERNAL BETA-CHANNEL V2 READINESS GATE

- Data UTC: 2026-03-04T07:35:01.734468+00:00
- Readiness: **READY**
- Verdict: **TRUE_EXTERNAL_BETA_CHANNEL_V2_READY_STRICT**
- pass_count: 8/8

## Counts
- n_pairs_total: 4000
- n_events_total: 12
- n_pre/n_post: 2706/1294
- interventions_used: 5

## Feature Diagnostics
- frac_f_std_eq0: 0.0003
- frac_f_autoc1_eq0: 0.0003
- frac_f_switch_eq0: 0.0568
- frac_f_slope_eq0: 0.0003
- median_abs_autoc1: 0.3129
- median_f_std: 0.1134
- median_local_neigh_n: 8.00

## Flags
- externality_ok: True
- manifest_roles_ok: True
- schema_pairs_ok: True
- schema_events_ok: True
- regime_split_ok: True
- power_floor_ok: True
- feature_non_degenerate_ok: True
- intervention_diversity_ok: True

## Required Next Step
- RUN_V2_BLIND_INTERVENTION_AND_TRIAD_VALIDATION

## Artifacts
- JSON: `report_qw2015_true_external_beta_channel_v2_readiness_gate.json`
