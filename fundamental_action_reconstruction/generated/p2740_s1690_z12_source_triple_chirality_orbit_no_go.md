# P2740/S1690 Z12 source-triple chirality orbit no-go

Status: `P2740_Z12_SOURCE_TRIPLE_CHIRALITY_ORBIT_NO_GO`

## Finite orbit audit
- ordered_distinct_triples=1320
- unordered_triples=220
- positive_ordered_triples=660
- negative_ordered_triples=660
- translation_unordered_orbit_count=19
- affine_unordered_orbit_count=9
- translation_orbits_with_nonzero_signed_sum=0
- affine_orbits_with_nonzero_signed_sum=0

## Theorem statement
The Z12 ordered-triple chirality is nonzero on every ordered distinct triple and is balanced globally (660 positive, 660 negative).  After forgetting source labels to translation or affine unordered orbits, every orbit has signed sum zero, because each unordered triple contains three positive and three negative orderings.  Thus the object supplies a real chiral sign only after choosing an ordered source triple; current artifacts provide no strict source localizer or P2721 polarity coupling for that choice.

## Recommendation
Do not promote ordered-triple chirality by itself: it is a real nonzero pointwise sign, but without a strict source-localizer for one ordered triple and an exported P2721 polarity-coupling theorem it remains label/order premise data.  The next proof-grade move must either construct that localizer-and-coupling theorem for this exact triple-chirality object, or pivot to a different strict signed observable with nonzero orbit-safe signed value; otherwise preserve the P2697-P2740 no-new-live-frontier certificate.
