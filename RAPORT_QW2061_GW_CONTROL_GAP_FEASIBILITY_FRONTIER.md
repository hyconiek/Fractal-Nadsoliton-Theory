# RAPORT QW-2061: GW CONTROL-GAP FEASIBILITY FRONTIER

- Data UTC: 2026-03-04T02:08:03.322422+00:00
- Verdict: **GW_CONTROL_GAP_FEASIBLE_IN_CURRENT_LINEAR_FEATURE_SPACE**
- n_candidates: 1771
- pass_count_all: 255
- soft_count_no_gap: 597

## Locked Flavor Context (from QW-2060)
- CKM/PMNS mean rel%: 11.867/9.386

## Best Row (primary score)
- weights: {'w_max_abs_corr': 0.0, 'w_mean_abs_corr': 0.65, 'w_corr_at_0ms': 0.2, 'w_corr_at_10ms': 0.15000000000000002}
- GW auc/adv/sep/gap: 0.8196/0.3300/0.002052/0.001277
- all_pass: True
- soft_pass_no_gap: True

## Best Global Control-Gap Row
- weights: {'w_max_abs_corr': 0.0, 'w_mean_abs_corr': 0.0, 'w_corr_at_0ms': 0.6000000000000001, 'w_corr_at_10ms': 0.4}
- GW auc/adv/sep/gap: 0.5136/0.2233/0.000092/0.000007

## Best Control-Gap Under Soft Constraints
- weights: {'w_max_abs_corr': 0.05, 'w_mean_abs_corr': 0.55, 'w_corr_at_0ms': 0.35000000000000003, 'w_corr_at_10ms': 0.05}
- GW auc/adv/sep/gap: 0.7607/0.3577/0.002172/0.001261

## Required Next Step
- LOCK_BEST_FEASIBLE_GW_WEIGHTS_WITH_LOCKED_FLAVOR_CONTEXT

## Artifacts
- JSON: `report_qw2061_gw_control_gap_feasibility_frontier.json`
