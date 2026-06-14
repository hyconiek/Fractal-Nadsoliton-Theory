# P2727/S1677 orientation-transition law equivariance and polarity no-go

Status: `P2727_ORIENTATION_TRANSITION_LAW_EQUIVARIANCE_POLARITY_NO_GO`

## Finite law enumeration
- checked_transition_rows=1152
- law_count=4
- equivariant_law_names=['preserve_orientation', 'flip_orientation']
- equivariant_polarity_selecting_law_count=0
- premise_polarity_selecting_law_names=['collapse_to_plus_orientation', 'collapse_to_minus_orientation']
The only inversion-equivariant source-independent laws are preserve and flip.  Preserve has zero jump; flip has balanced +4/-4 jumps.  The laws that select one nonzero polarity collapse to a chosen orientation and fail inversion equivariance, so they are premise selectors.

## Recommendation
Do not continue source-independent orientation-law variants.  A next admissible proof-grade move must introduce a source-dependent but non-premise strict invariant that breaks inversion equivariance with a computable signed value and an exported P2721 polarity-coupling theorem; otherwise preserve the P2697-P2727 no-new-live-frontier certificate.
