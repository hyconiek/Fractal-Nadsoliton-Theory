# P3047/S1997 memory-lag polarity source-law obstruction

Status: `P3047_MEMORY_LAG_POLARITY_SOURCE_LAW_OBSTRUCTION_BOUNDED_NO_EXPORT`

## Finite certificate
- source domain rows: `5`
- trivial source rows: `4`
- trivial-source equivariant maps: `0`
- odd-source equivariant maps: `2`
- candidate law rows: `3`
- accepted polarity source-law rows: `0`
- satisfied proof obligations: `3/6`
- strict polarity source law exported: `False`

## Decision
P3047 proves the finite source-domain obstruction for memory-lag polarity.  Aut-trivial source data have zero equivariant maps to the inversion-odd memory-lag sign torsor.  An inversion-odd source domain has the right representation type and two equivariant maps, but current artifacts export no nonzero strict odd source value, so no unique polarity-selection law is available.

## Recommendation
Do not replay positive sign, lag-2 winner, or even commutator magnitudes as a polarity source.  A next admissible move must supply one concrete nonzero strict inversion-odd source value coupled to the memory-lag sign, or pivot to a different genuinely new typed object; otherwise preserve the P3044-P3047 bounded no-export certificate.
