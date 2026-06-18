# P2842/S1792 exchangeable edge-pair measure localization candidate audit

Status: `P2842_EXCHANGEABLE_EDGE_PAIR_MEASURE_LOCALIZATION_CANDIDATE_NO_GO_NO_COUPLING_NO_CLOSURE`

## Candidate
The S_16-exchangeable probability measure on unordered vertex pairs with mass split by edge/non-edge occupancy on the fixed 16-node 4-regular carrier.

## Finite carrier check
- decoded_graph_count=16828
- distinct_edge_counts=[32]
- measure_constant_on_full_carrier=True

## Premise result
- accepted_as_strict_localization_or_coupling_object=False
- missing_premises=['nonconstant_on_full_carrier', 'separates_or_refines_combined_functional', 'field_or_spacetime_support_exported', 'coupling_or_variational_chain_rule_exported']

## Acceptance
- accepted_as_finite_label_safe_measure_candidate=True
- accepted_as_strict_localization_or_coupling_object=False
- accepted_as_bounded_candidate_no_go=True

## Boundary
P2842 introduces exactly one new candidate after P2841: the S_16-exchangeable edge-pair probability measure on the fixed 16-node 4-regular carrier.  It is label-gauge invariant and finite, with edge density 32/120=4/15.  However every decoded graph has the same edge count 32, so the measure is constant on the full 16,828-graph carrier; it does not refine the combined witness and exports no field/spacetime support or coupling/variational chain rule.  It is accepted only as a finite label-safe measure candidate, not as a strict localization/coupling object.

## Recommendation
Do not replay exchangeable edge-count measures or other carrier-constant summaries.  A next admissible move needs a genuinely nonconstant, label-safe strict localization object with field support and a coupling/chain-rule theorem; otherwise preserve the P2841-P2842 no-new-live-frontier/no-localization boundary.
