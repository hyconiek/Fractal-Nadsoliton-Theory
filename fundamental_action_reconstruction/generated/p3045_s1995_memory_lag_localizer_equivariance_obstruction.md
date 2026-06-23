# P3045/S1995 memory-lag localizer equivariance obstruction

Status: `P3045_MEMORY_LAG_LOCALIZER_EQUIVARIANCE_OBSTRUCTION_BOUNDED_NO_EXPORT`

## Finite certificate
- oriented lag rows: `11`
- translation-invariant rows: `11`
- inversion-flip rows: `11`
- localizer candidate rows: `4`
- accepted chart-independent lag-localizer rows: `0`
- Aut-compatible candidate rows: `0`
- satisfied proof obligations: `4/7`
- chart-independent lag localizer exported: `False`

## Decision
P3045 constructs the oriented lag torsor for the P3044 memory-lag commutator and verifies that integrated scores are cyclic-origin invariant.  The same computation shows the obstruction: Aut inversion sends each oriented lag score to the opposite lag with opposite sign, so receiver winners are chart-polarity choices rather than a nonpremise lag/localizer theorem.

## Recommendation
Do not promote translation-invariant lag scores or lag-2 winners to selector closure.  A next proof-grade move may attack exactly one remaining P3044 premise: either a strict nadsoliton source law for the memory-lag commutator/lag polarity, or an explicit coupling theorem from this signed lag torsor to a selector torsor/readout row.
