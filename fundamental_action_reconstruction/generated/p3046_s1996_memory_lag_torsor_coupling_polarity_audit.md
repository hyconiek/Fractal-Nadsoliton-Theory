# P3046/S1996 memory-lag torsor coupling-polarity audit

Status: `P3046_MEMORY_LAG_TORSOR_COUPLING_POLARITY_AUDIT_BOUNDED_NO_EXPORT`

## Finite certificate
- sign torsor size: `2`
- orientation torsor size: `2`
- candidate coupling rows: `2`
- Aut-equivariant coupling rows: `2`
- polarity-selected rows: `0`
- accepted selector/readout coupling rows: `0`
- source/readout rows: `3`
- accepted source/readout rows: `0`
- satisfied proof obligations: `3/6`
- selector/readout coupling exported: `False`

## Decision
P3046 finds real representation-level progress: the signed memory-lag torsor has exactly two Aut-equivariant maps to the orientation/selector torsor.  But the maps are opposite coupling-polarity choices, and current artifacts export no strict law selecting one polarity or installing it as a nonpremise selector/readout row.

## Recommendation
Do not promote the two conditional coupling maps to selector closure.  A next proof-grade move must supply a strict polarity-selection/source law for the memory-lag sign, or pivot to a different new typed object; otherwise preserve the P3044-P3046 bounded no-export certificate.
