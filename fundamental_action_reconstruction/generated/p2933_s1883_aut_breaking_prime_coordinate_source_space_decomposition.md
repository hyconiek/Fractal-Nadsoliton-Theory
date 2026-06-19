# P2933/S1883 Aut-breaking prime-coordinate source-space decomposition

Status: `P2933_AUT_BREAKING_PRIME_COORDINATE_SOURCE_SPACE_DECOMPOSITION_NO_SOURCE_LAW`

## Decomposition certificate
- product equations: `30`
- product rank/nullity: `6` / `5`
- Aut equations: `44`
- combined rank/nullity: `11` / `0`
- prime-coordinate basis count: `5`
- all basis vectors break Aut invariance: `True`
- accepted candidates: `0`

## Boundary
P2933 decomposes the exact additive value-source space after P2932.  Product additivity leaves five prime-coordinate degrees of freedom, while Aut(Z12)-invariance kills all five.  Therefore any nonzero L_p source must be an explicit strict Aut-breaking value law; current coordinate choices, external logs, and residue labels are not such a law.

## Recommendation
Construct a concrete Strict_AutBreaking_PrimeCoordinate_Source_Law that derives one nonzero vector in the five-dimensional prime-coordinate space from strict nadsoliton data and then run finite additivity, symmetry-breaking provenance, delta/eta, and beta/eta coupling tests.  If no such law is supplied, preserve the P2929-P2933 no-new-live-frontier boundary.
