# P3023/S1973 kernel-dissipation time-order candidate obstruction

Status: `P3023_KERNEL_DISSIPATION_TIME_ORDER_CANDIDATE_OBSTRUCTION_NO_CLOSURE`

## Finite certificate
- strict descent chain edges / total: `11/11`
- cyclic reset is strict descent: `False`
- U(12)-equivariant rows / total: `1/4`
- accepted as strict time-order with physical unit: `False`

## Decision
A new typed time-order candidate was constructed from K_strict_gate dissipation: K(d) strictly decreases along the finite chain 1->2->...->12.  The obstruction is that this chain is chart-dependent, has a cyclic reset at 12->1, is not U(12)-equivariant, and has no strict physical unit theorem for the directed step.

## Recommendation
Do not promote kernel-dissipation label chains to physical time.  A next proof-grade move may attack exactly one missing theorem for this object: a strict chart/selector source for the integer order or an independent physical tick theorem; otherwise pivot to a different genuinely new typed object while preserving the P3017-P3023 no-closure boundary.
