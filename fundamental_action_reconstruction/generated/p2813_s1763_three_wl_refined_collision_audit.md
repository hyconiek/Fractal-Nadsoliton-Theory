# P2813/S1763 three-WL refined collision audit

Status: `P2813_THREE_WL_HISTOGRAM_RESIDUAL_COLLISION_OBSTRUCTION_NO_CLOSURE`

## Candidate invariant
stable 3-dimensional Weisfeiler-Lehman ordered-triple color histogram, applied inside P2812 residual collision classes

## Exact quotient counts
- decoded_graph_count=16828
- p2812_two_wl_refined_class_count=16749
- p2812_collision_graph_count=158
- three_wl_computed_graph_count=158
- three_wl_refined_class_count=16749
- three_wl_collision_class_count=79
- three_wl_collision_graph_count=158
- three_wl_max_class_size=2
- remaining_defect_canonical_minus_three_wl=79
- defect_reduction_vs_p2812=0

## Decision
- accepted_as_three_wl_residual_obstruction_audit=True
- accepted_as_complete_canonical_source_carrier=False
- accepted_as_strict_source_law_or_coupling=False

## Boundary
P2813 applies exactly one stronger typed invariant, stable 3-WL ordered-triple color histograms, only to the 158 P2812 residual-collision graphs.  It is a finite negative obstruction for this histogram-level 3-WL candidate: the quotient remains at 16,749 classes, the defect remains 79, and all 79 collision pairs survive.  Therefore 3-WL histogram data do not export a complete canonical source carrier, strict spectral source law, or typed coupling theorem to K or L_total.

## Recommendation
Do not replay histogram-level 3-WL as a source carrier.  The next honest move should add exactly one different typed ingredient targeted at the 79 surviving pairs, preferably an exact automorphism/orbit partition or canonical-label-derived orbit signature, and only after actual separation run a separate source-law/coupling acceptance audit.  Keep role-bearing L_total, bridge closure, role transfer, selector closure, and ToE blocked.
