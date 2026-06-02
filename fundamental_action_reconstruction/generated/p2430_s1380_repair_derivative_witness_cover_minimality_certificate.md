# P2430 S1380: repair derivative witness-cover minimality certificate

Status: `PASS_REPAIR_DERIVATIVE_WITNESS_COVER_MINIMALITY_NO_GATE_DISCHARGE_NO_CLOSURE`

## Finite facts

- Global minimal covers: `[['source_obligation_discharge', 'chi11_source_export', 'qw2191_selector_discharge', 'role_transfer_audit_license', 'role_bearing_ltotal_export']]`.
- Minimal cover sizes by target: `{'bridge_source_ready': 1, 'selector_source_ready': 2, 'role_transfer_ready': 1, 'role_bearing_ltotal_ready': 1, 'toe_ready': 5}`.
- Global covering rows: `1`.
- Global proper failures: `31`.

## Hard limits

- Witness-cover minimality identifies theorem targets, not theorem-gate discharge.
- The unique ToE/global cover still contains all five missing theorem gates.
- Current APD value evidence and chi11 cut mechanism cover no theorem-gate witness requirement.
- No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'five_gates': True, 'five_targets': True, 'unique_global_all_five': True, 'unique_toe_all_five': True, 'selector_pair': True, 'global_one_cover': True, 'global_31_failures': True, 'toe_one_cover': True, 'toe_31_failures': True, 'uncovered_distribution_binomial': True, 'current_value_evidence_covers_no_theorem_gate': True, 'p2429_witness_rows_inherited': True, 'p2429_pairs_inherited': True, 'p2429_toe_inherited': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
