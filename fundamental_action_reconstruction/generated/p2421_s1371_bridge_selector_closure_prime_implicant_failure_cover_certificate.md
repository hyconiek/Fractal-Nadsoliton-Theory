# P2421 S1371: bridge-selector closure prime-implicant/failure-cover certificate

Status: `PASS_UNIQUE_PRIME_IMPLICANT_FAILURE_COVER_NO_GATE_DISCHARGE`

## Finite facts

- Closure gates: `7`.
- Assignments: `128`.
- True masks: `1`.
- ANF terms / degree: `1` / `7`.
- Prime implicant masks: `[127]`.
- Failure-cover literals: `7`.
- Current missing gates: `['source_obligation_discharge', 'chi11_source_export', 'qw2191_selector_discharge', 'role_transfer_audit_license', 'role_bearing_ltotal_export']`.

## Hard limits

- The unique prime implicant is a closure obligation, not a proof that its gates are discharged.
- The failure cover blocks every proper subset from ToE promotion.
- No source, chi11, QW-2191, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'seven_gates': True, 'all_128_assignments': True, 'one_true_mask': True, 'full_mask_127': True, 'single_anf_term': True, 'anf_degree_seven': True, 'unique_prime_implicant': True, 'seven_cnf_units': True, 'seven_failure_literals': True, 'all_127_proper_masks_fail': True, 'seven_nearest_misses': True, 'all_derivatives_decisive': True, 'current_missing_five_gates': True, 'current_repair_distance_five': True, 'p2420_minimal_mask_inherited': True, 'p2420_repair_distance_inherited': True, 'p2420_subcube_failure_inherited': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
