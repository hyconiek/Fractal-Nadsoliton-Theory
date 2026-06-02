# P2431 S1381: admissible next theorem-target antichain certificate

Status: `PASS_ADMISSIBLE_NEXT_TARGET_ANTICHAIN_NO_GATE_DISCHARGE_NO_CLOSURE`

## Finite facts

- Candidate rows size <= 2: `15`.
- Admissible candidates: `6`.
- Admissible singletons: `[['source_obligation_discharge'], ['chi11_source_export'], ['qw2191_selector_discharge']]`.
- Minimal readiness-complete antichain: `[['source_obligation_discharge'], ['chi11_source_export', 'qw2191_selector_discharge']]`.

## Hard limits

- The antichain is a target-selection fork, not theorem-gate discharge.
- Role-transfer and role-bearing L_total remain inadmissible first moves.
- The selector fork still requires both chi11 source export and QW-2191 discharge.
- No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'candidate_rows_15': True, 'admissible_candidates_6': True, 'inadmissible_candidates_9': True, 'admissible_singletons_expected': True, 'inadmissible_singletons_expected': True, 'minimal_readiness_antichain_expected': True, 'max_completed_coverage_antichain_expected': True, 'p2430_global_cover_inherited': True, 'p2430_selector_pair_inherited': True, 'role_transfer_not_first': True, 'ltotal_not_first': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
