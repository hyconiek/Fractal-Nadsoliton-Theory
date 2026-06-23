# P3038/S1988 viscous-retarded-chiral selector candidate obstruction

Status: `P3038_VISCOUS_RETARDED_CHIRAL_SELECTOR_CANDIDATE_OBSTRUCTION_NO_SELECTOR_EXPORT`

## Finite certificate
- operator rows: `4`
- finite nonzero operator rows: `4`
- accepted operator source rows: `0`
- P3037 features present: `6`
- plus branch score: `-0.007800396949201`
- minus branch score: `0.007800396949201`
- branch separation abs: `0.015600793898403`
- satisfied proof obligations: `4/8`
- strict selector mechanism exported: `False`

## Decision
The integrated candidate operator is real and finite: it combines memory/viscosity, c-retardation anisotropy, an inversion-odd chiral projection, and a unit/readout coupling slot, yielding nonzero opposite branch scores.  This is progress beyond hint collection, but the sign/localizer, path anisotropy, physical unit coupling, and nonproxy variational/readout theorem remain unsourced.  Therefore no strict selector mechanism is exported.

## Recommendation
Do not promote the branch-separating score alone to selector closure. The next proof-grade move should attack exactly one missing source premise of this integrated operator: either a nonpremise chiral sign/localizer theorem, a sourced retardation path-anisotropy theorem, or a physical unit/readout coupling theorem.
