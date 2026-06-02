# P2429 S1379: repair derivative nearest-miss witness table certificate

Status: `PASS_REPAIR_DERIVATIVE_WITNESS_TABLE_NO_GATE_DISCHARGE_NO_CLOSURE`

## Finite facts

- Witness rows: `69`.
- Witness count by target: `{'bridge_source_ready': 16, 'role_bearing_ltotal_ready': 16, 'role_transfer_ready': 16, 'selector_source_ready': 16, 'toe_ready': 5}`.
- Target/gate pairs with witnesses: `10`.
- ToE nearest-miss rows: `5`.

## Hard limits

- Derivative witness rows certify essential blockers, not theorem-gate discharge.
- The five ToE nearest misses remain five missing theorem obligations.
- The selector witness rows still require both chi11 source export and QW-2191 discharge.
- No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'truth_table_32': True, 'witness_rows_69': True, 'witness_counts_by_target': True, 'target_gate_pairs_10': True, 'toe_nearest_miss_5': True, 'toe_one_per_gate': True, 'toe_missing_before_one': True, 'selector_only_pair': True, 'ltotal_only_ltotal_export': True, 'p2428_derivatives_inherited': True, 'p2428_toe_inherited': True, 'p2428_selector_inherited': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
