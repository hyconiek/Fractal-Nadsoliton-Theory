# P3015/S1965 unit-action-compatible time observable functor obstruction

Status: `P3015_UNIT_ACTION_COMPATIBLE_TIME_OBSERVABLE_FUNCTOR_OBSTRUCTION_NO_CLOSURE`

## Finite certificate
- orbit count: `6`
- unit compatibility rows/failures: `48/0`
- successor rows / bad successor orbits: `6/4`
- accepted as time observable generator: `False`

## Decision
A genuine U(12)-unit-compatible, observer-independent strict-kernel observable functor was constructed by orbit-averaging.  The bounded obstruction is that the orbit quotient destroys the clock successor: same-orbit labels can advance to different successor orbits, so the functor is not a directed time observable generator and cannot install EOM/Hamiltonian or ToE closure.

## Recommendation
Do not replay orbit-average quotient observables as time-arrow closure.  The next proof-grade move should attack exactly one missing successor/evolution atom: construct a strict unit-compatible directed successor/semigroup on the quotient, or pivot to a unit-bearing action/EOM source for one typed observable; keep selector, bridge, role-transfer, L_total, observed-physics, and ToE closure blocked until that single atom passes finite acceptance.
