# P3037/S1987 selector-mechanism hint sheaf acceptance matrix

Status: `P3037_SELECTOR_MECHANISM_HINT_SHEAF_ACCEPTANCE_MATRIX_NO_SELECTOR_EXPORT`

## Finite certificate
- hint source rows: `4`
- required features: `6`
- features with some coverage: `6`
- glue profiles: `15`
- feature-complete profiles: `1`
- accepted selector-mechanism profiles: `0`
- exported mechanism rows: `0`
- strict selector mechanism exported: `False`

## Decision
The repository contains real selector hints outside a single familiar schema: viscosity/damping/memory, c-retardation with anisotropic anchor dependence, inversion-odd/chiral branch-source needs, and external unit/readout-coupling needs. The finite hint sheaf covers all six required feature atoms in some glue profiles, but zero source rows export a non-premise complete selector mechanism; therefore feature coverage is not selector closure.

## Recommendation
Do not force the selector into a familiar prechosen schema, but also do not promote hints to closure. The next proof-grade move should construct one new integrated candidate operator combining memory/viscosity, retardation anisotropy, an inversion-odd signed value, and explicit unit/readout coupling, then run this sheaf acceptance test on that single candidate.
