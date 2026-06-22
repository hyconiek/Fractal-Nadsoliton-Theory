# P3014/S1964 time-readout observable-generator candidate obstruction

Status: `P3014_TIME_READOUT_OBSERVABLE_GENERATOR_CANDIDATE_BOUNDED_NO_GO`

## Finite certificate
- explicit input/output types: `True`
- observer-independent formula: `True`
- unit equivariance rows/failures: `48/24`
- arrow rows passed/total: `3/4`
- accepted as time observable generator: `False`

## Decision
A concrete observer-independent time-readout formula was constructed from the strict kernel finite difference, but the finite test blocks it as a ToE-grade time observable: it is not U(12)-unit compatible, does not install an EOM/Hamiltonian, and cannot be promoted to a directed physical time arrow without the missing selector/unit/action atoms.

## Recommendation
Do not replay finite-difference kernel-clock observables.  The next proof-grade move should attack one different atom: either construct a unit-action-compatible observable functor for a single readout row, or provide a genuine unit-bearing action/EOM source for the already explicit time observable; keep selector, bridge, role-transfer, L_total, and ToE closure blocked unless that atom passes a finite acceptance test.
