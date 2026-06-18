# P2847/S1797 target-independent volume-form/unit-source audit

Status: `P2847_TARGET_INDEPENDENT_VOLUME_FORM_UNIT_SOURCE_AUDIT_NO_CLOSURE`

## Carrier check
- decoded_graph_count=16828
- coverage_ok=True

## Candidate rows
### uniform_vertex_density
- distinct_raw_total_mass_count=1
- raw_total_mass_range=[16, 16]
- distinct_normalized_probability_profile_count=1
- accepted_as_target_independent_unit_volume_source=False
- missing_premises=['nonconstant_density_available', 'unit_dimension_source', 'canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'locality_covariance', 'coupling_coefficient_rule', 'variational_chain_rule']
### triangle_plus_one_density
- distinct_raw_total_mass_count=1
- raw_total_mass_range=[16, 16]
- distinct_normalized_probability_profile_count=1
- accepted_as_target_independent_unit_volume_source=False
- missing_premises=['nonconstant_density_available', 'unit_dimension_source', 'canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'locality_covariance', 'coupling_coefficient_rule', 'variational_chain_rule']
### four_cycle_plus_one_density
- distinct_raw_total_mass_count=36
- raw_total_mass_range=[36, 232]
- distinct_normalized_probability_profile_count=2866
- accepted_as_target_independent_unit_volume_source=False
- missing_premises=['target_independent_total_volume', 'unit_dimension_source', 'canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'locality_covariance', 'coupling_coefficient_rule', 'variational_chain_rule']
### triangle_square_plus_one_density
- distinct_raw_total_mass_count=36
- raw_total_mass_range=[36, 232]
- distinct_normalized_probability_profile_count=2866
- accepted_as_target_independent_unit_volume_source=False
- missing_premises=['target_independent_total_volume', 'unit_dimension_source', 'canonical_vertex_to_field_support', 'spacetime_pullback_formula', 'locality_covariance', 'coupling_coefficient_rule', 'variational_chain_rule']

## Boundary
P2847 finds that finite probability normalization is available for the tested label-safe densities, and nonconstant normalized profiles exist.  However, variable raw masses make nonconstant density-weighted total volumes carrier-dependent, while the uniform count is only a dimensionless counting measure.  No candidate exports a unit dimension source, canonical vertex-to-field support, spacetime pullback, coupling coefficient rule, or variational chain rule.

## Recommendation
Do not replay finite probability normalization as a unit-bearing volume form.  The next admissible move should isolate the alternative P2845 premise: one target-independent coupling coefficient/unit source for a single named density, with an explicit source law and units.  If that source law is not supplied, pivot to exactly one concrete kernel bridge atom with exported source premise, or preserve the no-new-live-frontier certificate.
