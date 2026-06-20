# P2963/S1913 typed mediator/functor no-go

Status: `P2963_TYPED_MEDIATOR_FUNCTOR_NO_GO_NO_STRICT_EXPORT`

## Mediator certificate
- mediator count: `8`
- equalizing mediators: `['max_nonzero']`
- accepted typed mediators: `[]`
- strict typed mediator exported: `False`
- acceptance matrix rows/accepted: `32/1`

## Lay summary
P2963 tests the typed-mediator escape route explicitly.  It finds that max_nonzero can equalize K and C, but only by erasing the support/provenance mismatch identified in P2962.
No strict mediator is exported: every mediator that preserves the mismatch keeps K and C unequal, while every equalizer is too coarse to be a typed provenance functor.

## Recommendation
Do not replay typed-mediator scalar collapses, max-coordinate equalization, coefficient-quotient exchange, target_sum cuts, primitive equal-summand, K+C decompositions, or beta-scale normalization.  A next proof-grade move must either construct an actual unit-bearing nonproxy coupling that receives the P2938 aggregate without needing K/C artifact exchange, or introduce a new typed structural object richer than the current K/C split; otherwise pivot outside the ratio-package lane while preserving the P2929-P2963 no-strict-export boundary.
