# P2817/S1767 structural observable source-law obstruction audit

Status: `P2817_STRUCTURAL_OBSERVABLE_SOURCE_CANDIDATE_OBSTRUCTED_NO_COUPLING_NO_CLOSURE`

## Candidate
Q_struct(G)=(edge density, degree histogram, distance histogram, exact 4-cycle count)

## Normalization candidate
edge density is normalized by C(16,2); histograms are finite count measures; 4-cycle count is dimensionless

## Typed map candidate
G -> Q_struct(G) -> dimensionless structural coefficient vector for K or L_total

## Finite counts
- decoded_graph_count=16828
- q_struct_class_count=228
- q_struct_collision_class_count=165
- q_struct_collision_graph_count=16765
- q_struct_max_class_size=788
- remaining_defect_canonical_minus_q_struct=16600

## Acceptance
- accepted_as_structural_observable_obstruction_audit=True
- accepted_as_source_law_coupling=False
- accepted_as_bounded_candidate_rejection=True

## Boundary
P2817 tests exactly one non-ordinal normalized structural observable candidate after P2816: Q_struct=(edge density, degree histogram, distance histogram, exact 4-cycle count).  The finite computation covers all 16,828 graphs, but Q_struct collapses the complete P2815 carrier into far fewer classes and leaves residual collisions.  It also exports no strict graph-source law or typed K/L_total coupling theorem.  Therefore P2817 is a bounded structural-observable obstruction, not dynamics or closure.

## Recommendation
Do not replay ordinal ranks or this low-order structural observable as source-law evidence.  The next honest move is exactly one richer non-ordinal formula with an actual typed coupling ansatz, preferably a local edge-toggle variational response energy whose value is computed directly from the full P2815 separating response data and accompanied by a falsifiable graph-to-K or graph-to-L_total normalization/coupling theorem.  If no such formula is supplied, preserve the no-coupling boundary.
