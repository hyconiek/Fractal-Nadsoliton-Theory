# P2811/S1761 local motif-refined source candidate audit

Status: `P2811_LOCAL_MOTIF_REFINED_SOURCE_CANDIDATE_AUDIT_REFINES_BUT_STILL_OBSTRUCTED_NO_CLOSURE`

## Candidate class
F(G)=Phi(chi_A_G, chi_A_complement_G, local_motif_distance_profile_G)

## Exact quotient counts
- decoded_graph_count=16828
- spectral_pair_class_count_from_p2810=16211
- refined_class_count=16691
- canonical_target_class_count_from_p2810=16828
- refined_collision_class_count=132
- refined_collision_graph_count=269
- refined_max_class_size=3
- remaining_defect_canonical_minus_refined=137
- defect_reduction_vs_spectral_pair_only=480

## Decision
- accepted_as_local_motif_refinement_audit=True
- accepted_as_complete_canonical_source_carrier=False
- accepted_as_strict_source_law_or_coupling=False

## Boundary
P2811 adds one richer local motif/distance profile beyond the P2810 spectral-pair data.  The refinement improves the quotient from 16,211 to 16,691 classes and reduces the defect from 617 to 137, but 132 collision classes covering 269 graphs remain.  Therefore this richer invariant is useful evidence but still not a complete canonical source carrier, strict spectral source law, or K/L_total coupling theorem.

## Recommendation
Attack the remaining 132 refined collision classes directly with one additional typed invariant, preferably an automorphism/orbit-aware motif invariant or a canonical-label-derived local orbit signature, and report whether it separates the 269 remaining graphs.  Do not promote local motif refinement to K/L_total or ToE while any refined collisions remain and no typed variational coupling is exported.
