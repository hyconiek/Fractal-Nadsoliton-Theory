# NOC-004 — Parent-versus-reference ledger

- Actual parents: two separate **unchanged** q12 labelled LOO HB laws with A7 and diagnostic A7+.01(P3-P4). Original denominator N and slot rate 1 unchanged.
- Auxiliary reflected laws: two model-specific local Gibbs laws, restricting only total-count exits from the inherited domains; never called physical barriers.
- Start: actual iid-uniform labels at model time zero, then one shared deterministic polynomial observation origin `t_N` late enough for both local mixing windows. Prior outside-domain and genuine exits are kept.
- Transfer: couple actual and reflected *entire path* on `[t_N,t_N+M_N]`; do not transport unbounded F/R moments by TV and do not postselect surviving trajectories.
- Bound: for each model j, original wrong-classification probability is ≤stationary-reflected wrong-event probability + `mu0(D^c)+mixing-TV+genuine-exit-probability`. Each error occurs once.
- Final decision: uses only `M_N` observed block aggregates, no future `Y_M` endpoint, no full n, no soft coordinate, no outcome-dependent block.
- Limitation: the effective upper bound is asymptotic and conditional on prior MT domain/mixing/trapping proofs; neither finite-size N0, runtime calibration, error prefactor, physical sample requirement nor exponent optimality is claimed.
