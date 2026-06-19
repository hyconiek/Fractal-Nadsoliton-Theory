# P2952/S1902 torsion-character aggregate weight-source obstruction

Status: `P2952_TORSION_CHARACTER_AGGREGATE_WEIGHT_SOURCE_OBSTRUCTION_NO_STRICT_PROVENANCE`

## Weight-source certificate
- weight domain: `[1, 2, 3, 4]`
- family rows: `16`
- all rows product-additive by construction: `True`
- target weight pair count in domain: `1`
- target weight pairs: `[[1, 1]]`
- target vector forces equal weights: `True`
- strict equal-weight source theorem exported: `False`
- P2951 torsion-character provenance atom discharged: `False`

## Boundary
P2952 attacks the P2951 torsion-character provenance atom by splitting the P2938 aggregate into kernel-excess and character-negativity ingredients.  The exact P2938 vector forces the equal-weight pair a=b=1 inside V_p(a,b)=a*K_p+b*C_p, but current artifacts do not export a strict nadsoliton theorem selecting that equal-weight aggregation over the positive weight family.  The provenance atom therefore remains undischarged.

## Recommendation
Do not add more bounded weight scans or alternative linear aggregate weights.  A next proof-grade move must export an actual strict equal-weight/source theorem for K_p+C_p, or attack a different P2951 atom with a concrete theorem object; otherwise pivot outside the ratio-package lane and preserve the P2929-P2952 no-strict-provenance boundary.
