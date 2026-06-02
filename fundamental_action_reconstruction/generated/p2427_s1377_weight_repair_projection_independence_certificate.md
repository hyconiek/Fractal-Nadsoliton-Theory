# P2427 S1377: weight-repair projection independence certificate

Status: `PASS_WEIGHT_REPAIR_PROJECTION_INDEPENDENCE_NO_GATE_DISCHARGE_NO_CLOSURE`

## Finite facts

- Weight assignments: `144`.
- Repair assignments: `32`.
- Product assignments: `4608`.
- Repair readiness true counts: `{'bridge_source_ready': 16, 'selector_source_ready': 8, 'role_transfer_ready': 16, 'role_bearing_ltotal_ready': 16, 'toe_ready': 1}`.
- Missing-count distribution: `{'0': 1, '1': 5, '2': 10, '3': 10, '4': 5, '5': 1}`.
- Tables factor: `True`.

## Hard limits

- Projection independence is not a theorem-gate discharge.
- Weight winner labels have empty support on repair readiness predicates.
- Repair subsets have empty support on weighted frontier winner predicates.
- No source, selector, QW-2191, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'weight_grid_144': True, 'repair_subcube_32': True, 'product_4608': True, 'weight_counts_match_p2425': True, 'repair_readiness_counts_expected': True, 'missing_distribution_binomial': True, 'gate_inclusions_half_subcube': True, 'all_tables_factor': True, 'conditional_readiness_identical': True, 'conditional_missing_identical': True, 'empty_weight_support_for_readiness': True, 'empty_repair_support_for_winner': True, 'p2426_product_inherited': True, 'p2426_winner_counts_inherited': True, 'p2426_toe_count_inherited': True, 'p2426_no_gate_reduction_inherited': True, 'weighted_choice_discharges_no_gate': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
