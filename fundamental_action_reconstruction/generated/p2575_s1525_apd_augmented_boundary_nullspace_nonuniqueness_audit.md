# P2575/S1525 APD augmented-boundary nullspace nonuniqueness audit

Status: `P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Augmented basis: `['V=prod_{d=0}^{11}(x-d)', 'x*V', 'x^2*V']`.
- Boundary matrix rank/nullity: `2/1`.
- Targets audited: `3`.
- All targets solved by augmented basis: `True`.
- Nullspace changes off-node values while preserving boundaries: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

Adding `V`, `x*V`, and `x^2*V` repairs the P2574 two-endpoint compatibility obstruction, but the resulting boundary matrix has rank `2` and nullity `1`.  Therefore the augmented family can satisfy the audited endpoint slopes while still carrying a null direction that preserves all finite APD nodes and boundary slopes but changes off-node values.  Solvability by extra vanishing modes is not a source theorem.

## Recommended next honest step

Do not repair two-endpoint APD incompatibility by freely adding extra vanishing modes and then claim a source theorem. The next honest step is to derive the admissible APD function space and boundary data from the strict action; otherwise the boundary problem can be made solvable but remains nonunique through a nullspace.

## Negative controls

No APD augmented boundary source, vanishing-mode source, boundary nullspace selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`e2c29b968bb68290beeda8c7f9d20f3497e7ec53b80eeeebc9acfa35d2492278`
