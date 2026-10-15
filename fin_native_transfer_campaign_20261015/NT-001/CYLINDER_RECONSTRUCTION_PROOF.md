# NT-001 — Full native temporal reflection positivity

**Disposition:** `EXACT_FINITE_REVERSIBLE_CHAIN_THEOREM` for the stationary *reflected local count chain*. The same finite stationary identity applies separately to the global Gibbs count chain, but the two stationarities must not be conflated. This is not a new microscopic process.

## 1. Parent and reversible measure

Let `E_N=D_N` be the admitted connected local count graph and `n=(n_0,...,n_11)` with `sum n_i=N`. For `i!=j`, `n_i>0`, and `n'=n-e_i+e_j in D_N`, the row generator is

`L_D(n,n') = n_i softmax_j[(g_N/N) A (n-e_i)]`, `L_D(n,n)=-sum_{n'!=n} L_D(n,n')`.

The same formula without domain suppression defines the unreflected count process. Self-refresh attempts give zero generator increment. Under the already accepted constant diagonal of A7, the local Gibbs count law is

`pi_D(n) ∝ (N!/(12^N prod_i n_i!)) exp[(g_N/(2N)) n^T A n] 1_{n∈D_N}`.

For `m=n-e_i`, the forward/reverse softmax denominators agree because the reverse departure background is also `m`. Thus

`L_D(n,n')/L_D(n',n) = [n_i/(n_j+1)] exp[(g_N/N)((A m)_j-(A m)_i)]`.

The Gibbs weight ratio has the identical value: the energy difference equals `(g_N/N)((Am)_j-(Am)_i)` since `A_jj=A_ii`. Hence `pi_D(n)L_D(n,n')=pi_D(n')L_D(n',n)` for all retained edges. The chain is finite, conservative and reversible, with strictly positive invariant pi_D on D_N. Unreflected global Gibbs has its **own** equilibrium and does not concentrate in D_N.

## 2. Reflection-positive two-sided cylinders

Let X_t be the two-sided stationary Markov process with invariant pi_D, transition `P_t=e^{t L_D}` and reverse-time reflection `theta(t)=-t`. For any bounded complex functional `F` of finitely many future times `0<=t_1<...<t_k`, define `Theta F` by pulling those times to the past and conjugating coefficients. By the conditional independence of past and future given **the full X_0 count**, and by detailed balance (so the reversed past has the same conditional transition laws as the future),

`E_pi[conj(F o theta) F] = sum_n pi_D(n) |E[F | X_0=n]|^2 >=0`.

For a linear combination `F=sum_l c_l F_l`, the same identity holds with `E[F|X_0]` replaced by the corresponding linear combination. Equivalently, `G_lm=E[conj(Theta F_l)F_m]` is a positive semidefinite Gram matrix. For cylinders excluding time zero, approximation at times decreasing to zero is licensed by finite-state stochastic continuity. There is no assumption that the collective coordinate `z_N` alone is Markov.

## 3. Hilbert quotient and transfer

Define `K F=E[F|X_0]` in `L²(pi_D)`. The reflection inner product is `<F,G>_RP = <K F,K G>_pi`. Quotient by `ker K` and complete. Functions of X_0 embed isometrically and span all `L²(pi_D)` because the state space is finite, so the reconstructed Hilbert space is precisely this space (as an abstract Hilbert space). Positive future-time translation `tau_s` satisfies `K(tau_s F)=P_s(KF)` by the Markov property. Reversibility implies `P_s=P_s*`, and the Markov property plus spectrum(L_D) in `(-∞,0]` makes `P_s` a positive selfadjoint contraction. Its nonnegative generator is `-L_D`.

Let `D_pi=diag(pi_D)`; the counting-space transform `K_N= -D_pi^(1/2) L_D D_pi^(-1/2)` is a symmetric positive semidefinite matrix with ground vector `sqrt(pi_D)`. For `v=sqrt(pi_D) u`, the edge decomposition is

`<v,K_N v> = (1/2) sum_{n,n'} pi_D(n)L_D(n,n') |u(n')-u(n)|²`.

This is a **native transfer operator** made from actual interacting rates; it is not asserted to be conserved physical energy, a microscopic spacetime operator, or an independent quantum dynamics.

## 4. Initial-state and scope boundaries

* Stationary `pi_D` is a reflected analysis process, not the actual iid initialization.
* MT late-time labelled/path TV and departure control can transfer any **fixed bounded** path-cylinder expectation into the actual unreflected iid parent, with its error explicitly inherited. It does not produce an exact stationary full-space transfer operator for the iid parent.
* Full global Gibbs is reversible separately but its dominant equilibrium is not the critical symmetric local sector.
* No physical quantum Hilbert-space preparation or detector follows from an abstract L² construction.
