# P2428 S1378: repair-readiness ANF derivative certificate

Status: `PASS_REPAIR_READINESS_ANF_DERIVATIVE_NO_GATE_DISCHARGE_NO_CLOSURE`

## Finite facts

- Gate count: `5`.
- Truth-table size: `32`.
- ANF polynomials: `{'bridge_source_ready': 'source_obligation_discharge', 'selector_source_ready': 'chi11_source_export * qw2191_selector_discharge', 'role_transfer_ready': 'role_transfer_audit_license', 'role_bearing_ltotal_ready': 'role_bearing_ltotal_export', 'toe_ready': 'source_obligation_discharge * chi11_source_export * qw2191_selector_discharge * role_transfer_audit_license * role_bearing_ltotal_export'}`.
- Derivative edge counts: `{'bridge_source_ready': {'source_obligation_discharge': 16, 'chi11_source_export': 0, 'qw2191_selector_discharge': 0, 'role_transfer_audit_license': 0, 'role_bearing_ltotal_export': 0}, 'selector_source_ready': {'source_obligation_discharge': 0, 'chi11_source_export': 8, 'qw2191_selector_discharge': 8, 'role_transfer_audit_license': 0, 'role_bearing_ltotal_export': 0}, 'role_transfer_ready': {'source_obligation_discharge': 0, 'chi11_source_export': 0, 'qw2191_selector_discharge': 0, 'role_transfer_audit_license': 16, 'role_bearing_ltotal_export': 0}, 'role_bearing_ltotal_ready': {'source_obligation_discharge': 0, 'chi11_source_export': 0, 'qw2191_selector_discharge': 0, 'role_transfer_audit_license': 0, 'role_bearing_ltotal_export': 16}, 'toe_ready': {'source_obligation_discharge': 1, 'chi11_source_export': 1, 'qw2191_selector_discharge': 1, 'role_transfer_audit_license': 1, 'role_bearing_ltotal_export': 1}}`.

## Hard limits

- ANF/derivative structure is a blocker map, not a theorem-gate discharge.
- The selector monomial still requires both chi11 source export and QW-2191 discharge.
- The ToE monomial still requires all five missing theorem gates.
- No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'truth_table_32': True, 'five_gates': True, 'single_prime_each': True, 'selector_pair': True, 'toe_all_five': True, 'toe_derivative_nearest_misses': True, 'bridge_derivative_only_source': True, 'selector_derivatives_pair_only': True, 'ltotal_derivative_only_ltotal': True, 'p2427_product_inherited': True, 'p2427_factor_inherited': True, 'p2427_weight_support_empty_inherited': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
