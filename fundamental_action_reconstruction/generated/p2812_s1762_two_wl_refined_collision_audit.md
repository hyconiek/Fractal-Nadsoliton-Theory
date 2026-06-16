# P2812/S1762 two-WL refined collision audit

Status: `P2812_TWO_WL_REFINED_COLLISION_AUDIT_REFINES_BUT_STILL_OBSTRUCTED_NO_CLOSURE`

## Candidate invariant
stable 2-dimensional Weisfeiler-Lehman ordered-pair color histogram, applied inside P2811 residual collision classes

## Exact quotient counts
- decoded_graph_count=16828
- p2811_refined_class_count=16691
- p2811_refined_collision_graph_count=269
- two_wl_computed_graph_count=269
- two_wl_refined_class_count=16749
- two_wl_collision_class_count=79
- two_wl_collision_graph_count=158
- two_wl_max_class_size=2
- remaining_defect_canonical_minus_two_wl=79
- defect_reduction_vs_p2811=58

## Decision
- accepted_as_two_wl_refinement_audit=True
- accepted_as_complete_canonical_source_carrier=False
- accepted_as_strict_source_law_or_coupling=False

## Boundary
P2812 applies one additional typed invariant, stable 2-WL ordered-pair color histograms, only to the 269 P2811 residual-collision graphs.  It improves the refined quotient from 16,691 to 16,749 classes and reduces the remaining defect from 137 to 79, but 79 two-WL collision classes covering 158 graphs remain.  Thus 2-WL refinement is useful computational evidence, not a complete canonical source carrier, strict source law, or K/L_total coupling.

## Recommendation
Target the remaining 79 two-WL collision pairs with exactly one stronger invariant, preferably a 3-WL/k-WL signature, exact automorphism-orbit partition, or canonical-label-derived orbit signature.  Keep K/L_total, role transfer, bridge closure, selector closure, and ToE promotion blocked until separation and an explicit typed variational coupling are exported.
