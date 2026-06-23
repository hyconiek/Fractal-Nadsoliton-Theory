# P3039/S1989 chiral projector sign/localizer obstruction

Status: `P3039_CHIRAL_PROJECTOR_SIGN_LOCALIZER_OBSTRUCTION_NO_SOURCE_EXPORT`

## Finite certificate
- receiver rows: `4`
- finite nonzero rows: `4`
- inversion-odd rows: `1`
- accepted nonpremise localizer rows: `0`
- translation phase orbit size: `12`
- phase/polarity projector count: `12`
- satisfied proof obligations: `4/8`
- nonpremise chiral sign/localizer exported: `False`

## Decision
The P3038 chiral projector has a real finite inversion-odd torsor, and unit 11 flips its sign.  However translations move the phase origin, inversion exchanges polarity, and K/memory/absolute-K weighted receivers choose only chart-relative representatives.  Therefore no nonpremise chiral sign/localizer theorem is exported.

## Recommendation
Do not replay sine/chiral phase receivers as a selector source.  The next proof-grade move should attack one different P3038 missing source premise: either a sourced retardation path-anisotropy theorem for the c-retarded split, or a physical unit/readout coupling theorem for the integrated density.  Continue only with a genuinely new strict source law; otherwise preserve P3038-P3039 as a branch-separating-but-unsourced certificate.
