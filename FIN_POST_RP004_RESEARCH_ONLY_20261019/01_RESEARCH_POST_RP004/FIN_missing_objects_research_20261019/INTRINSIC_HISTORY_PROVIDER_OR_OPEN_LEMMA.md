# FIN: necessary history-provider object and what existing equations actually generate

Date: 2026-10-19. Scope: F01/F03/F06. No new source-law postulate adopted.

## A. A mathematically concrete history object constructible after a law is selected

On strict Z12, the supplied weights determine `A=sI-W`. Declare, as a **conditional choice of kinetics**, the continuous-time Markov generator `Q=-A`. `Q_ij=W_ij>0` for `i!=j`, `Q_ii=-s`, stationary distribution `pi_i=1/12`, and transition matrices `P_t=e^{tQ}`. For any `0<=t1<...<tn` the consistent finite-dimensional cylinder probabilities are

`H(i1,t1;...;in,tn)=pi_i1 product_{k=1}^{n-1} [P_{t_{k+1}-t_k}]_{i_k,i_{k+1}}`.

The finite-state generator theorem ensures a stationary path law with these cylinder probabilities. Thus **the history object exists constructively** on the chosen clock, but the step `Q=-A` is a kinetic choice, not a necessity inferred from static correlations. It cannot be used to claim that F01 is closed. `exp(-itA)` and `cos(t sqrt A)` remain mathematically distinct evolutions of the same operator.

The observed two-time joint matrix is `H_t=diag(pi)P_t`. A FULL, externally/independently supplied and temporally calibrated `H_t` for all small t yields

`Q = lim_{t->0+} [diag(pi)^(-1)H_t-I]/t`, if this limit exists.

For a fixed t>0, if additionally `P_t` is known to be reversible with strictly positive eigenvalues and a compatible reversible embedding is admitted, the real symmetric logarithm recovers its reversible generator: `Q=(1/t)log P_t`. This is **identification from independently provided dynamic data**, not source from kernel. Coarse readouts generally do not identify all entries.

## B. Exact pointwise quadratic-variation source in the jump class

Let the observable state of a finite 12-label jump process be the one-hot vector `X=e_i` when the label is `i`. At state i, a jump to j has increment `v_ij=e_j-e_i` and off-diagonal rate `q_ij`. The infinitesimal **increment second-moment matrix** (conditional, without prematurely averaging over pi) is

`Gamma_i = lim_{dt->0+} 1/dt E[(X(dt)-X(0))(X(dt)-X(0))^T | X(0)=e_i]`
`= sum_{j != i} q_ij v_ij v_ij^T`.

For every j != i, **`q_ij=(Gamma_i)_{jj}`**: in the j,j entry of a single-jump increment matrix only the i->j jump contributes. The diagonal is fixed by `q_ii=-sum_{j!=i}q_ij`. Hence the full pointwise Gamma_i family, **if independently supplied**, uniquely identifies Q in this labelled single-jump class. The formula is a transparent special-case reconstruction of the already accepted RC/AC inverse perspective, not a new source provider. An averaged quadratic variation or static pi does not supply Gamma_i.

A necessary compatibility certificate for a proposed independently supplied Gamma is:
  1. `Gamma_i= sum_j gamma_ij (e_j-e_i)(e_j-e_i)^T` with `gamma_ij>=0` and allowed support;
  2. `q_ij=gamma_ij` and `q_ii=-sum_j gamma_ij` define a conservative generator;
  3. If reversibility is claimed: `pi_i gamma_ij=pi_j gamma_ji`;
  4. If a temporal reading is claimed: one common clock maps its rates to the supplied history;
  5. For physical meaning: the state/readout and spacetime composition maps are separately proved, not inferred from Gamma by a change of notation.

## C. Why the currently available equilibrium does NOT yield the needed provider

For any fixed strictly positive 12-state stationary vector pi, all symmetric edge activities `c_ij=c_ji>0` define reversible rates `q_ij=c_ij pi_j`. Choosing `c_ij=2` reproduces the usual native N2 reset off-diagonal rates, but choosing `c_01=c_10=2(1+alpha)` with alpha>0 and all other c_ij=2 yields the **same** pi, symmetry under detailed balance, and different conditional Gamma and waiting times. It is not an overall rescaling of the clock. This is an independent elementary witness of the previously documented underdetermination, not a new universal no-go theorem.

Also the Markov path object constructed in A already assumes Q. Computing Gamma from that path and then claiming Gamma independently **generated Q** is a circle.

## D. The exact remaining missing object

The minimal meaningful fundamental input is a non-circular, typed `H_source` (or `Gamma_source`) with:
- its carrier, state and observation algebra specified independently of an already chosen Q or PDE;
- consistent histories / a law of intrinsic relational changes rather than merely an equal-time covariance;
- an independently justified comparison-event or time map, if model clock numbers are used;
- composition/coupling information and an explicit alternative generator that this object excludes;
- a derivation from already posited FIN assumptions that does **not** insert the answer as a hidden premise.

Such an object is **not present among the input equations of this work**. The present note pays exact algebraic construction and reconstructibility *conditional on a chosen law/provider*, not the existence of an independently justified fundamental provider. Therefore F01/F03/F06 remain OPEN_SOURCE.
