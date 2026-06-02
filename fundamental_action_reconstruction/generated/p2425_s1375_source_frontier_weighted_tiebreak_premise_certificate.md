# P2425 S1375: source-frontier weighted tie-break premise certificate

Status: `PASS_WEIGHTED_TIEBREAK_CONDITIONAL_NO_SOURCE_PREMISE_NO_DISCHARGE`

## Finite facts

- Weight grid assignments: `144`.
- Winner counts: `{'bridge_first_pareto': 108, 'bridge_selector_tie': 6, 'selector_pair_first_pareto': 30}`.
- Bridge-first condition: `w_selector < 2*w_bridge`.
- Selector-pair-first condition: `w_selector > 2*w_bridge`.
- Tie boundary: `w_selector = 2*w_bridge`.
- Dominated win count: `0`.

## Hard limits

- Weighted tie-break conditions are conditional on an extra source-cost premise.
- No internal weight/source-cost premise is exported by this certificate.
- No unique first source gate is selected in strict-core terms.
- No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'three_objective_classes': True, 'two_pareto_classes': True, 'grid_12_by_12': True, 'winner_counts_expected': True, 'tie_boundary_expected': True, 'bridge_condition_expected': True, 'selector_condition_expected': True, 'dominated_never_wins': True, 'mixed_dominated_symbolic': True, 'p2424_pareto_count_inherited': True, 'p2424_dominated_count_inherited': True, 'p2424_no_unique_inherited': True, 'no_internal_weight_premise': True, 'no_unique_first_gate': True, 'source_still_open': True, 'chi11_still_open': True, 'qw2191_still_open': True, 'role_transfer_still_blocked': True, 'ltotal_still_blocked': True, 'toe_still_open': True, 'fingerprint_stable': True}`
