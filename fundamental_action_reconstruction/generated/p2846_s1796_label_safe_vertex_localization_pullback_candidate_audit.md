# P2846/S1796 label-safe vertex localization/pullback candidate audit

Status: `P2846_LABEL_SAFE_VERTEX_LOCALIZATION_PULLBACK_CANDIDATE_NO_GO_NO_CLOSURE`

## Carrier check
- decoded_graph_count=16828
- coverage_ok=True

## Candidate rows
### uniform_vertex_measure
- distinct_profile_count=1
- singleton_profile_count=0
- max_profile_class_size=16828
- nonconstant_on_full_carrier=False
- accepted_as_strict_localization_pullback=False
- missing_premises=['nonconstant_on_full_carrier', 'canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'target_independent_units_or_volume_form', 'locality_covariance', 'variational_chain_rule', 'coupling_coefficient_rule']
### triangle_count_vertex_measure
- distinct_profile_count=1
- singleton_profile_count=0
- max_profile_class_size=16828
- nonconstant_on_full_carrier=False
- accepted_as_strict_localization_pullback=False
- missing_premises=['nonconstant_on_full_carrier', 'canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'target_independent_units_or_volume_form', 'locality_covariance', 'variational_chain_rule', 'coupling_coefficient_rule']
### four_cycle_count_vertex_measure
- distinct_profile_count=2877
- singleton_profile_count=1420
- max_profile_class_size=236
- nonconstant_on_full_carrier=True
- accepted_as_strict_localization_pullback=False
- missing_premises=['canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'target_independent_units_or_volume_form', 'locality_covariance', 'variational_chain_rule', 'coupling_coefficient_rule']
### triangle_square_joint_vertex_measure
- distinct_profile_count=2877
- singleton_profile_count=1420
- max_profile_class_size=236
- nonconstant_on_full_carrier=True
- accepted_as_strict_localization_pullback=False
- missing_premises=['canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'target_independent_units_or_volume_form', 'locality_covariance', 'variational_chain_rule', 'coupling_coefficient_rule']
### local_motif_wl_vertex_measure
- distinct_profile_count=672
- singleton_profile_count=395
- max_profile_class_size=11968
- nonconstant_on_full_carrier=True
- accepted_as_strict_localization_pullback=False
- missing_premises=['canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'target_independent_units_or_volume_form', 'locality_covariance', 'variational_chain_rule', 'coupling_coefficient_rule']

## Boundary
P2846 finds real nonconstant label-safe vertex-density profiles on the full 16,828 graph carrier, but every candidate remains an anonymous finite vertex measure.  None exports a canonical vertex-to-field support, spacetime pullback formula, target-independent volume/unit form, locality/covariance theorem, coupling coefficient rule, or variational chain rule.  Thus the P2845 localization premise remains blocked.

## Recommendation
Do not replay anonymous vertex-profile localization.  The next admissible move must add one genuinely new bridge from anonymous vertices to field support: either a strict coordinate/support theorem for x_G/rho_G(x), or a target-independent volume form tied to one named density.  If no such object is supplied, pivot to the alternative P2845 route: one target-independent coupling coefficient/unit source or one concrete kernel bridge atom with source premise.
