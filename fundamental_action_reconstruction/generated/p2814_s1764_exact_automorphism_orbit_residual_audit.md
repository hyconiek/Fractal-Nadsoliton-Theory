# P2814/S1764 exact automorphism-orbit residual audit

Status: `P2814_EXACT_AUTOMORPHISM_ORBIT_RESIDUAL_OBSTRUCTION_NO_CLOSURE`

## Candidate invariant
exact automorphism-group order and vertex-orbit partition signature, applied inside P2813 residual collision classes

## Exact quotient counts
- decoded_graph_count=16828
- p2813_refined_class_count=16749
- p2813_collision_graph_count=158
- automorphism_orbit_computed_graph_count=158
- automorphism_orbit_refined_class_count=16771
- automorphism_orbit_collision_class_count=57
- automorphism_orbit_collision_graph_count=114
- automorphism_orbit_max_class_size=2
- remaining_defect_canonical_minus_automorphism_orbit=57
- defect_reduction_vs_p2813=22
- automorphism_group_order_histogram_on_residual={1: 112, 2: 40, 4: 6}
- orbit_size_multiset_histogram_on_residual={'(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)': 112, '(2, 2, 2, 2, 2, 2, 2, 2)': 40, '(4, 4, 4, 4)': 6}
- truncated_automorphism_search_count=0

## Decision
- accepted_as_exact_automorphism_orbit_obstruction_audit=True
- accepted_as_complete_canonical_source_carrier=False
- accepted_as_strict_source_law_or_coupling=False

## Boundary
P2814 applies exactly one different typed ingredient, exact automorphism-group/order and vertex-orbit partition signatures, only to the 158 P2813 residual-collision graphs.  The computation is exact and untruncated; it improves the quotient from 16,749 to 16,771 classes and reduces the defect from 79 to 57.  However 57 collision pairs covering 114 graphs remain, so automorphism/orbit data still do not export a complete canonical source carrier, strict spectral source law, or typed coupling theorem to K or L_total.

## Recommendation
Do not replay automorphism/orbit signatures as closure.  The next honest move should add exactly one non-automorphism typed ingredient targeted at the 57 surviving pairs, preferably a small-subgraph extension/deletion response profile or a typed graph-to-K/L_total coupling ansatz with a falsifiable acceptance matrix.  Keep role-bearing L_total, bridge closure, role transfer, selector closure, and ToE blocked.
