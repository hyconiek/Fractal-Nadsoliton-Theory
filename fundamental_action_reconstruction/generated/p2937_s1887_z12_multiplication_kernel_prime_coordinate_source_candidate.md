# P2937/S1887 Z12 multiplication-kernel prime-coordinate source candidate

Status: `P2937_Z12_MULTIPLICATION_KERNEL_PRIME_COORDINATE_SOURCE_CANDIDATE_REJECTED_PARTIAL`

## Candidate certificate
- vector [L2,L3,L5,L7,L11]: `[1, 2, 0, 0, 0]`
- nonzero prime coordinates: `2`
- zero prime coordinates: `[5, 7, 11]`
- product pairs d*e<=11: `29`
- product additivity defects: `0`
- accepted strict prime-log source: `False`

## Boundary
P2937 constructs an explicit finite Z12 multiplication-kernel candidate instead of replaying coordinate scans.  The candidate is product-additive and intrinsic to multiplication on Z/12Z, but it only sources nonzero coordinates for 2 and 3; unit primes 5, 7, and 11 receive zero.  It also lacks strict nadsoliton provenance and delta/eta plus beta/eta coupling, so it is rejected as a strict L_p source.

## Recommendation
Do not continue multiplication-kernel variants unless a new strict theorem couples unit-prime information to nonzero 5/7/11 coordinates and to delta/beta damping.  The next admissible move is a genuinely new unit-sensitive strict source object, or preservation of the P2929-P2937 no-new-live-frontier boundary.
