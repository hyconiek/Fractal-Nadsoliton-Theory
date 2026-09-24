# OP-005 protocol theorem / counterexample note

## Sequential ledger theorem

Let stages `k=1,...,K` be reached according to an adapted protocol. Let `E_k` denote a false promoted claim at stage `k`. If for every null parameter and every history with positive null probability,

`P(E_k | F_(k-1), reach k) <= delta_k`,

then

`P(union_k E_k) <= sum_k delta_k`.

Proof: write `P(E_k and no earlier false claim)` as an expectation of its conditional probability on the history/reach event and bound it by `delta_k P(reach k and no earlier false claim) <= delta_k`; sum over the first false stage. The argument permits adaptive stopping and predictable branch selection.

## E-value composition condition

If nonnegative factors `e_k` satisfy `E[e_k|F_(k-1)]<=1`, then `M_n=product_(k<=n)e_k` is a nonnegative supermartingale and is safe under optional stopping in the usual theorem scope.

Marginal expectations alone are insufficient. Counterexample: `e_1=e_2=2B` for `B~Bernoulli(1/2)`. Then `E[e_1]=E[e_2]=1` but `E[e_1 e_2]=2>1`.

## Required data contract

The current repository does not prove a shared-data joint law covering calibration selection, rank testing and channel classification simultaneously. Therefore the default admissible implementation uses disjoint scoring blocks after an independent calibration/training block. Reuse is allowed only if a future theorem establishes simultaneous confidence/e-process validity for that reuse.
