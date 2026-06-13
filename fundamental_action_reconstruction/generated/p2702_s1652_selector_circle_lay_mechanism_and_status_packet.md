# P2702/S1652 selector-on-circle lay mechanism and status

Status: `P2702_SELECTOR_CIRCLE_LAY_MECHANISM_AND_STATUS_PACKET_NO_FALSE_PASS`

## Lay explanation
- **What is the selector?** A selector is a rule that chooses one representative from many symmetric possibilities.  On a 12-node circle it is like choosing one arrow around the circle instead of treating all equivalent arrows as interchangeable.  Status: `mechanism_real_as_concept`
- **How does it break symmetry and pick one direction?** A premise selector can mark +1 as the arrow.  In the toy model this reduces the residual unit symmetry from [1,5,7,11] to [1], so one direction is chosen.  Status: `works_if_the_selector_is_given_as_a_premise`
- **What about convexity on the circle?** The 'convexity break on a circle' is best read as cutting a cyclic/flat degeneracy into one branch or sector.  It is intuitive branch selection, not a strict convexity theorem exported here.  Status: `lay_metaphor_or_branch_cut_not_a_new_strict_convexity_theorem`
- **Why does current repo not claim final closure?** Because P2699/P2700 show Aut-invariant internal information does not choose +1 over -1, and P2701 finds no new strict-sourced provider.  The earlier direction support is real but premise-based.  Status: `no_false_promotion`

## Next honest step
The next proof-grade move cannot be another explanation or Aut-invariant replay.  It must construct an actual new strict-sourced symmetry-breaking provider with non-premise provenance, or introduce a different new typed object outside closed lanes.  Without that, keep the P2697-P2702 no-new-live-frontier certificate.
