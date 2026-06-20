# P2962/S1912 K/C exchange-symmetry source obstruction

Status: `P2962_KC_EXCHANGE_SYMMETRY_SOURCE_OBSTRUCTION_NO_STRICT_EXPORT`

## Obstruction certificate
- coefficient quotient candidate available: `True`
- support cardinality match: `False`
- nonzero multiset match: `False`
- typed provenance equivalence exported: `False`
- strict artifact symmetry source exported: `False`
- unit-bearing nonproxy coupling exported: `False`
- strict K/C exchange source exported: `False`
- acceptance matrix rows/accepted: `64/1`

## Lay summary
P2962 sharpens P2961: the quotient exchange is mathematically clean, but the actual P2938 K and C artifacts are not currently typed-isomorphic.
No strict exchange source yet.  The obstruction is concrete: support sizes 2 vs 3, nonzero multisets [1,2] vs [2,2,2], and distinct provenance labels block an artifact-level K/C symmetry on current data.

## Recommendation
Do not replay coefficient-quotient exchange, target_sum=9 cuts, primitive equal-summand convention, K+C decompositions, or beta-scale normalization.  A next proof-grade move must introduce a new typed mediator/functor that makes K and C comparable without erasing their support/provenance mismatch, or construct a unit-bearing nonproxy coupling that does not require artifact-level K/C exchange symmetry; otherwise pivot outside the ratio-package lane while preserving the P2929-P2962 no-strict-export boundary.
