# P2749/S1699 minimal inversion-odd source coupling-polarity audit

Status: `P2749_MINIMAL_ODD_SOURCE_COUPLING_POLARITY_GAP`

## Finite coupling audit
- orientation_reversing_units=[7, 11]
- all_set_maps_count=4
- equivariant_map_count=2
- equivariant_maps=[{'-1': -1, '1': 1}, {'-1': 1, '1': -1}]
- unique_coupled_map_count_after_p2721_polarity=2

## Theorem statement
A minimal inversion-odd signed source is the correct representation type requested by P2748: there are exactly 2 Aut(Z12)-equivariant maps from the odd source sign to the orientation torsor.  But those two maps are opposite coupling polarities, and composing with P2721 polarity exchanges them.  Therefore supplying an abstract odd sign is not enough: a further strict law must choose the source sign and the coupling polarity.

## Recommendation
Do not promote an abstract minimal inversion-odd source to selector closure.  P2749 shows the representation type is admissible but leaves exactly two opposite equivariant couplings, exchanged by P2721 polarity.  The next proof-grade move must provide a concrete strict source sign value and a theorem selecting one coupling polarity; otherwise pivot to a different typed object or preserve the P2697-P2749 no-new-live-frontier certificate.
