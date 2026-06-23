# P3040/S1990 retarded path-anisotropy source obstruction

Status: `P3040_RETARDED_PATH_ANISOTROPY_SOURCE_OBSTRUCTION_NO_SOURCE_EXPORT`

## Finite certificate
- receiver rows: `4`
- finite nonzero rows: `4`
- accepted path-anisotropy source rows: `0`
- rho grid rows: `12`
- best grid rows: `2`
- path sign flip verified: `True`
- satisfied proof obligations: `4/8`
- sourced retardation path anisotropy exported: `False`

## Decision
The P3038 retarded split is finite and nonzero, and the path-sign torsor is explicit: reversing the sign sends Delta to -Delta.  But current receivers only insert rho/sign, tune them to the branch readout, or use chart-order gradient/curvature diagnostics.  No strict path geometry exports a positive rho, parallel/perpendicular retarded sectors, or a chart-independent path localizer.

## Recommendation
Do not replay rho-grid tuning, cyclic gradients, memory curvature, or inserted parallel/perpendicular labels as a sourced retarded geometry.  The remaining direct P3038 premise is a physical unit/readout coupling theorem for the integrated density; otherwise preserve P3038-P3040 as a finite branch-separating but unsourced selector-candidate certificate.
