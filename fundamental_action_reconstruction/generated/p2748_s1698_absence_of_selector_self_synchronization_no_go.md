# P2748/S1698 absence-of-selector self-synchronization no-go

Status: `P2748_ABSENCE_OF_SELECTOR_SELF_SYNCHRONIZATION_NO_GO`

## Finite synchronization audit
- orientation_reversing_units=[7, 11]
- singleton_absence_map_count=2
- singleton_absence_equivariant_map_count=0
- bit_absence_map_count=4
- bit_absence_equivariant_map_count=0
- odd_signed_source_equivariant_map_count=2

## Theorem statement
An Aut(Z12)-invariant absence/no-selector datum has trivial action.  Because units 7 and 11 reverse the orientation torsor, every map from a fixed absence state (or a trivial absence bit) to the selector torsor fails equivariance: the finite counts are 0 equivariant maps from the singleton and 0 from the trivial bit.  Equivariant synchronization becomes possible only after replacing absence by a new inversion-odd signed source; that has 2 equivariant maps, but it is exactly the missing new sign, not information about absence alone.

## Recommendation
Do not promote 'information about absence of selector' to a selector.  P2748 proves that an Aut-invariant absence/no-selector datum has zero equivariant maps to the orientation torsor; synchronization requires an added inversion-odd signed source, which is precisely the missing new object.  The next proof-grade move must either construct that concrete inversion-odd signed source with a P2721 coupling theorem, or pivot to a different typed object; otherwise preserve the P2697-P2748 no-new-live-frontier certificate.
