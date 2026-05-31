# Legacy-linear vs strict-nonlinear damping compression separation certificate

Status: `legacy-linear-torsion-denominator-separated-from-strict-nonlinear-d-eta-compression-no-bridge-theorem`

## Separation summary

- `domain`: `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]`
- `eta`: `1.8`
- `eta_minus_one`: `0.8`
- `legacy_beta_tors`: `0.01`
- `strict_beta`: `1.0`
- `required_gamma_values_strictly_increase`: `True`
- `no_single_linear_gamma_matches_two_distinct_positive_nodes`: `True`
- `legacy_beta_tors_matches_no_positive_strict_node`: `True`
- `best_l2_linear_gamma`: `5.562387511264106`
- `best_l2_residual_l2_norm`: `26.820335473174516`
- `best_l2_residual_max_abs`: `13.718051778840163`
- `legacy_beta_tors_residual_min`: `0.99`
- `legacy_beta_tors_residual_max`: `74.79431440274533`
- `required_gamma_at_d1`: `1.0`
- `required_gamma_at_d11`: `6.809483127522302`
- `required_gamma_spread_d11_minus_d1`: `5.809483127522302`
- `d11_strict_over_legacy_denominator_ratio`: `68.38226522769848`
- `component_gap_records_compression_missing`: `True`
- `guardrail_records_legacy_incomplete`: `True`
- `necessity_marks_damping_shape_critical`: `True`
- `exact_damping_monotone_inherited`: `True`
- `compression_ontology_identifies_missing_characteristic`: `False`
- `t15_names_damping_compression_obligation`: `True`
- `full_bridge_theorem_exported`: `False`

## Cross-checks

- `source_gap_and_guardrail_loaded`: `True`
- `required_gamma_strictly_increases`: `True`
- `no_single_linear_gamma_matches_two_nodes`: `True`
- `legacy_beta_tors_matches_no_positive_node`: `True`
- `best_linear_fit_still_has_residual`: `True`
- `damping_shape_critical_and_monotone_inherited`: `True`
- `bridge_obligation_named_but_not_closed`: `True`

All cross-checks pass: `True`

## Proof certificate

- `nonduplication_step`: rg was used to separate this linear-vs-nonlinear denominator no-go from prior damping monotonicity, inverse-branch, ontology, and component-gap reports.
- `algebraic_step`: If 1+gamma*d=1+d^(9/5) for d>0, then gamma=d^(4/5); because d^(4/5) is strictly increasing on positive d, no constant gamma can match two distinct positive nodes.
- `legacy_beta_tors_step`: The legacy beta_tors=0.01 member matches no positive strict node on d=1..11; the required gamma already equals 1 at d=1 and increases to 11^(4/5).
- `least_squares_step`: The best finite L2 linear denominator fit has gamma=5.562387511264, but its residual norm is 26.820335473175, so even the best constant linear torsion replacement is not exact on the audited domain.
- `bridge_meaning_step`: The strict d^eta denominator is therefore a genuine strict-side compression addition that must be supplied by the completion map; it is not hidden inside beta_tors*d.
- `theoretical_limit`: This is a finite/algebraic separation certificate for the damping-compression row, not a strict-source derivation of eta, beta, or the full legacy->strict bridge.

## Hard limits

- No raw identity K_legacy_ont == K_strict_gate is claimed.
- No derivation of eta=9/5, beta=1, or strict d^eta compression from nadsoliton dynamics is exported here.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer is licensed.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
