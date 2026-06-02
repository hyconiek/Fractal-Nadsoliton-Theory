# P2422 S1372: current missing-gate repair subcube certificate

Status: `PASS_CURRENT_REPAIR_SUBCUBE_PARTIAL_UNLOCKS_NO_GATE_DISCHARGE`

## Finite facts

- Missing gates: `5`.
- Repair rows: `32`.
- Proper repair failures: `31`.
- ToE repair count: `1`.
- Selector singleton unlock count: `0`.
- Singleton non-ToE unlock count: `3`.
- ToE-ready added gates: `['source_obligation_discharge', 'chi11_source_export', 'qw2191_selector_discharge', 'role_transfer_audit_license', 'role_bearing_ltotal_export']`.

## Hard limits

- Repair-subcube readiness is conditional on adding theorem gates; it does not discharge them.
- Selector-source readiness has no singleton repair and still requires chi11 source plus QW-2191 discharge.
- All 31 proper repair subsets remain ToE failures.
- No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'five_missing_gates': True, 'all_32_repair_rows': True, 'proper_31_failures': True, 'single_toe_repair': True, 'toe_requires_all_five': True, 'bridge_source_singleton': True, 'selector_pair_required': True, 'role_transfer_singleton': True, 'ltotal_singleton': True, 'toe_minimal_all_five': True, 'no_selector_singleton': True, 'three_singleton_non_toe_unlocks': True, 'p2421_missing_gates_inherited': True, 'p2421_repair_distance_inherited': True, 'p2421_prime_implicant_inherited': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
