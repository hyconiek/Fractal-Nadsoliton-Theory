# Legacy-to-strict damping parameter identifiability certificate

Status: `strict-beta-eta-finitely-identifiable-from-denominator-samples__source-and-legacy-map-open`

## Result

The strict denominator samples identify `beta=1` and `eta=9/5` on the
audited positive domain, and satisfy the fifth-power cover identity
`(S(d)-1)^5=d^9`.  This is not a source theorem.

## Summary

- `positive_domain`: `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]`
- `beta_identified_from_d1`: `True`
- `identified_beta`: `1`
- `identified_eta`: `9/5`
- `all_beta_recoveries_equal_1_given_eta`: `True`
- `max_beta_recovery_residual`: `0.0`
- `all_eta_log_recoveries_equal_9_over_5_given_beta`: `True`
- `max_eta_log_recovery_residual`: `4.440892098500626e-16`
- `exact_fifth_power_cover_identity_recorded`: `True`
- `candidate_grid_p_le_30_q_le_10_matching_count`: `1`
- `candidate_grid_unique_match`: `True`
- `required_linear_gamma_strictly_increases`: `True`
- `legacy_beta_tors_not_equal_any_required_gamma`: `True`
- `damping_separation_inherited`: `True`
- `exact_damping_calculus_inherited`: `True`
- `finite_diagonal_completion_map_inherited`: `True`
- `component_gap_damping_source_still_open`: `True`
- `strict_beta_eta_source_exported`: `False`
- `legacy_beta_tors_to_beta_eta_theorem_exported`: `False`
- `full_bridge_theorem_exported`: `False`

## Cross-checks

- `guardrail_records_strict_compression_missing_from_legacy`: `True`
- `beta_eta_recovered_from_samples`: `True`
- `finite_cover_identity_and_candidate_unique`: `True`
- `legacy_linear_not_parameter_match`: `True`
- `prior_damping_and_diagonal_reports_inherited`: `True`
- `no_source_or_full_bridge_claim`: `True`

## Proof certificate

- `grep_step`: rg was used to distinguish this beta/eta identifiability certificate from prior damping monotonicity and linear-vs-nonlinear separation reports.
- `beta_step`: At d=1, S(1)-1=beta*1^eta=beta, and the strict sample has S(1)-1=1, so beta is identified as 1 within the accepted strict denominator model.
- `eta_step`: For every d>=2, eta=log(S(d)-1)/log(d)=log(d^(9/5))/log(d)=9/5 once beta=1 is fixed.
- `cover_step`: Equivalently, eta=9/5 places the strict denominator increment on the finite algebraic cover (S(d)-1)^5=d^9.
- `candidate_step`: The rational grid p/q with p<=30 and q<=10 has exactly one matching candidate, 9/5, when beta is fixed by the d=1 sample and d=2 is checked.
- `legacy_linear_step`: The required linear gamma at node d is d^(4/5), which strictly increases, so no single beta_tors-style linear damping parameter supplies the strict denominator samples.
- `theoretical_limit`: This identifies beta/eta from accepted strict samples only; it is not a strict dynamical source for beta/eta, not beta_tors->chi_11, not legacy role transfer, and not ToE closure.

## Hard limits

- No strict dynamical derivation of beta=1 or eta=9/5 is exported.
- No beta_tors -> beta/eta theorem is claimed.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer is licensed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
