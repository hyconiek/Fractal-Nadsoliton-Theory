# MT-002 — High-probability long residence from the unchanged count generator

**Status:** analytic theorem conditional on MT-001's local domain/positive boundary. Applies to q=12 N∈12ℕ, |λ|≤Λ, the *actual* unit-attempt leave-one-out heat-bath chain, and the two frozen diagnostic MET/BAR rules at their own rates. No reflected process is used in the trapping calculation.

Let μ_{0,N}(n)=N!/(12^N ∏_i n_i!), the iid-uniform initial occupation distribution (not a stationary law). Define the UNNORMALIZED finite count measure

 W_N(n)=μ_{0,N}(n) exp[(g_N/(2N)) n^T A n].

For constant-diagonal A, the self-interaction correction in the underlying pair Gibbs energy is the same for every count type and cancels in normalized weights. Every frozen rate rule satisfies finite-N detailed balance with W_N. The primary LOO heat bath follows from conditional Gibbs single-copy refresh; MET and BAR from the edge energy difference exp((g_N/N)(A(n−e_i))_j−... ) with fixed prefactors.

Since A=XX^T⪰0 and g_N≥0, W_N(n)≥μ0_N(n) POINTWISE. As finite positive measures and using W_N P_t=W_N,

 μ0_N P_t ≤ W_N P_t=W_N          for every t≥0. (1)

This is invariant-*weight* domination, not an assertion that iid uniform is stationary. The finite count chain is conservative; applying a positive Markov kernel preserves order of finite measures.

For every type p=n/N the method-of-types upper bound yields

 μ0_N(n) ≤exp(−N D(p||u)), W_N(n)≤exp(−N F_{g_N}(p)). (2)

At most (N+1)^12 different count types exist. Let B_N={n∈D_N: an actual one-count edge from n to D_N^c exists}. The actual outgoing total count rate satisfies ∑_{n'≠n}L_N(n,n')≤2N: HB≤N, MET≤(11/12)N, BAR≤(22/12)N. The compensator for counts of exiting jumps gives

 P_iid(τ_D≤T)≤μ0_N(D_N^c)+ E_iid ∫_0^T 1_{n(t)∈B_N} ∑_{n'∉D_N}L_N(n(t),n')dt
 ≤μ0_N(D_N^c)+2NT∑_{n∈B_N}W_N(n)
 ≤μ0_N(D_N^c)+2NT(N+1)^12 exp(−Nb). (3)

The pre-jump state is in B_N. Equation (1) is used at **each** time t and no equilibrium initial law is inserted. The definition τ_D=0 on initial outside states accounts for the first term. There is no hidden postselection, no ignored invalid event and no old B4 barrier.

For fixed c with 0<c<b, set T_N=exp(cN). By (3), the second term tends to zero (polynomial factor times exp[−N(b−c)]). Along with MT-001 iid coverage,

 P_iid(τ_D>exp(cN))→1. (4)

This is a one-sided asymptotic residence guarantee. It is NOT a matching mean-exit asymptotic, a large-deviation quasipotential equality, a physical number of seconds, a frozen parity sign, an autonomous bit or a global phase-switching theorem. Because the example h and b are deliberately very small, no *practical* N lifetime is inferred.
