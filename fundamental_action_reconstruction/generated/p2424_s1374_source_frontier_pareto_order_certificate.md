# P2424 S1374: source-frontier Pareto order certificate

Status: `PASS_SOURCE_FRONTIER_PARETO_ORDER_NO_UNIQUE_SOURCE_GATE_NO_DISCHARGE`

## Finite facts

- Source-frontier gates: `3`.
- Admissible orders: `6`.
- Pareto orders: `4`.
- Dominated orders: `2`.
- Pareto vectors: `{'[1, 3]': 2, '[3, 2]': 2}`.
- Pareto classes: `{'bridge_first_pareto': 2, 'selector_pair_first_pareto': 2}`.

## Hard limits

- Pareto optimality ranks admissible proof orders but does not discharge a source theorem.
- Two source-frontier Pareto classes remain incomparable without an extra cost/source premise.
- No unique first source gate is selected internally.
- No role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'three_source_frontier_gates': True, 'six_admissible_orders': True, 'four_pareto_orders': True, 'two_dominated_orders': True, 'class_counts_expected': True, 'pareto_class_counts_expected': True, 'objective_vectors_expected': True, 'pareto_vectors_expected': True, 'no_unique_first_gate': True, 'p2423_admissible_orders_inherited': True, 'p2423_role_step_inherited': True, 'p2423_ltotal_step_inherited': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
