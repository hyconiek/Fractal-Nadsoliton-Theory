# RAPORT QW-1947: COUPLING COMPLETENESS + MASS OPERATOR FRONTIER

- Data UTC: 2026-03-03T06:53:57.092400+00:00
- Verdict: **COUPLING_AUDIT_AND_MASS_OPERATOR_FRONTIER_STRICT_FAIL**

## Strict (Baseline Assignment)
- strict_pass_exists: False
- best strict operator: O1_hard_linear_local14 | mean/max rel%=78.152/99.890 | complexity_k=1

## Exploratory (Assignment Pool from QW-1943)
- exploratory_pass_exists: False
- best exploratory operator: O1_hard_linear_local14 @ audit_best_mass_assignment | mean/max rel%=71.635/99.752

## Coupling Completeness Audit
- legacy QW-1939 feature coverage: 1/8 (amplitude only)
- tested union coverage: 8/8
- best strict operator feature coverage: 1/8
- never tested canonical features: NONE
- best strict missing canonical features: ['local_curvature', 'local_gradient', 'memory_cumulative', 'nonlocal_scale_mix', 'parity_structure', 'phase_offset', 'sign_oscillation']

## Pareto Frontier (Strict mean error vs complexity)
- O1_hard_linear_local14: mean/max rel%=78.152/99.890, k=1

## Required Next Step
- MASS_LINK_REWORK_REQUIRED_BEYOND_TESTED_DETERMINISTIC_OPERATOR_CLASSES

## Artifacts
- JSON: `report_qw1947_coupling_completeness_and_mass_operator_frontier.json`
