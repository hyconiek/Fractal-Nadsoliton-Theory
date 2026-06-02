# P2420 S1370: bridge-selector nonclosure reason matrix certificate

Status: `PASS_BRIDGE_SELECTOR_MECHANISM_NONCLOSURE_REASON_MATRIX_NO_TOE_CLOSURE`

## Direct answer

The current artifacts give a value-level APD bridge witness and locate a chi11 selector cut. They do not export the strict dynamic/source theorems, do not discharge QW-2191, do not perform the post-bridge role-transfer audit, and do not construct a role-bearing L_total; therefore ToE closure is blocked.

## Finite facts

- Closure gates: `7`.
- Closure assignments: `128`.
- Current true gates: `['apd_value_bridge_witness', 'chi11_phase_selector_cut_mechanism']`.
- Current mask: `5`.
- Minimal ToE masks: `[127]`.
- APD+selector-mechanism subcube failures: `31/32`.
- Current repair distance to ToE: `5`.

## Reasons

- R1_VALUE_BRIDGE_IS_NOT_SOURCE_BRIDGE
- R2_SELECTOR_MECHANISM_IS_NOT_SELECTOR_SOURCE
- R3_QW2191_IS_STILL_AN_INDEPENDENT_SELECTOR_OBSTRUCTION
- R4_ROLE_TRANSFER_AUDIT_IS_POST_BRIDGE_AND_NOT_AUTOMATIC
- R5_LTOTAL_EXPORT_IS_A_SEPARATE_ROLE_BEARING_STEP

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'seven_closure_gates': True, 'all_128_assignments': True, 'current_mask_is_apd_plus_mechanism': True, 'full_mask_127': True, 'unique_minimal_toe_mask': True, 'only_one_toe_true_mask': True, 'subcube_size_32': True, 'subcube_only_one_closure': True, 'subcube_31_failures': True, 'repair_distance_five': True, 'no_single_flip_closes_toe': True, 'five_reasons_recorded': True, 'apd_value_bridge_inherited': True, 'source_zero_inherited': True, 'chi11_cut_inherited': True, 'chi11_source_still_open': True, 'qw2191_still_open': True, 'bridge_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
