# P2426 S1376: weighted tie-break x repair-subcube nonpromotion product certificate

Status: `PASS_WEIGHTED_REPAIR_PRODUCT_NO_GATE_DISCHARGE_NO_CLOSURE`

## Finite facts

- Weight assignments: `144`.
- Repair assignments: `32`.
- Product assignments: `4608`.
- ToE-ready product rows: `144`.
- Proper repair failures: `4464`.
- Empty repair remaining-missing distribution: `{'5': 144}`.

## Hard limits

- Weight-side tie-breaks do not discharge theorem-gate repairs.
- Every proper repair subset remains a ToE failure for every positive weight assignment.
- The empty repair row keeps all five missing theorem gates under every weight winner.
- No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'weight_grid_144': True, 'repair_subcube_32': True, 'product_4608': True, 'toe_ready_144': True, 'proper_failures_4464': True, 'winner_counts_scaled': True, 'toe_counts_match_weight_grid': True, 'empty_repair_all_weights': True, 'empty_repair_all_five_missing': True, 'p2425_winner_counts_inherited': True, 'p2425_no_weight_premise_inherited': True, 'p2422_repair_count_inherited': True, 'p2422_toe_count_inherited': True, 'weighted_choice_does_not_reduce_missing_gates': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
