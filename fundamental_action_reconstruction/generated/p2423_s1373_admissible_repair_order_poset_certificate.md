# P2423 S1373: admissible repair-order poset certificate

Status: `PASS_ADMISSIBLE_REPAIR_ORDER_POSET_NO_GATE_DISCHARGE`

## Finite facts

- Missing gates: `5`.
- Precedence edges: `4`.
- Raw orders: `120`.
- Admissible orders: `6`.
- Rejected orders: `114`.
- Role-transfer step distribution: `{'4': 6}`.
- L_total step distribution: `{'5': 6}`.
- ToE step distribution: `{'5': 6}`.

## Hard limits

- Admissible order is not theorem discharge for any missing gate.
- Role-transfer remains after bridge/selector source gates and is not a shortcut around them.
- Role-bearing L_total remains last and cannot be promoted from earlier prefixes.
- No source, selector, role-transfer, L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'five_missing_gates': True, 'four_precedence_edges': True, 'all_120_raw_orders': True, 'six_admissible_orders': True, 'one_hundred_fourteen_rejected': True, 'role_transfer_only_step_four': True, 'ltotal_only_step_five': True, 'toe_only_step_five': True, 'selector_ready_step_distribution_expected': True, 'all_edges_necessary': True, 'p2422_repair_subcube_inherited': True, 'p2422_toe_repair_inherited': True, 'p2422_no_selector_singleton_inherited': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
