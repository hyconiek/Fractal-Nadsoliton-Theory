# RAPORT QW-2037: BETA DERIVATION BRANCH RESOLUTION AUDIT

- Data UTC: 2026-03-03T20:30:55.453526+00:00
- Readiness: **BETA_DERIVATION_GAP_OPEN**
- Verdict: **BETA_BRANCH_RESOLUTION_FAIL**
- pass_count: 2/5

## Target
- target beta (QW-2030): 1.170000

## Beta CI95
- opt beta q02.5/q50/q97.5: 0.550166 / 0.777247 / 0.990327
- resolved beta q02.5/q50/q97.5: 0.550166 / 0.751705 / 0.992875
- median rel objective increase (resolved): 0.000000

## Transfer Medians
- corr opt/resolved: 0.054142 / 0.054249
- gain opt/resolved: 0.000486 / 0.000504

## Flags
- target_beta_in_ci95_opt: False
- target_beta_in_ci95_resolved: False
- median_beta_distance_improved: False
- median_rel_obj_increase_le_0p25: True
- median_corr_not_worse: True

## Required Next Step
- EXTEND_MICRO_DYNAMICS_FOR_BETA_IDENTIFIABILITY

## Artifacts
- JSON: `report_qw2037_beta_derivation_branch_resolution_audit.json`
