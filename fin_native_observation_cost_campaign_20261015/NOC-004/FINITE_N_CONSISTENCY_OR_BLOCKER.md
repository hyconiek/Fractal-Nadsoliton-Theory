# NOC-004 — One consistent finite-parent bounded classifier, conditional theorem

For `N in 12N`, freeze **the exact unmodified** `k_N=floor(N^(15/16))`, `M_N=floor(N^(1/4))`, `alpha_N=k_N/N` and unit time grid. Eventually `1<=k_N<=N/2`; small-size examples where this fails are not used to alter the frozen asymptotic test. On the predetermined size-k block, measure only

`W_N(i)=k_N^-1/2 sum_{a in B_N} f3(x_a(t_N+i))`, `A_N=M_N^-1 sum_{i=0}^{M_N-1} cos(W_N(i))`.

For the two known *analysis predictions* `mu_{0,alpha_N}>mu_{1,alpha_N}`, choose the frozen midpoint `theta_N=(mu_0+mu_1)/2`. Declare `control` when `A_N<theta_N` and `A7` otherwise (ties go to A7). There is no input of the full counts or hidden state to this decision.

By NOC-002, `mu_0-mu_1>=c_*alpha_N` for sufficiently large N. Under either reflected stationary law, a wrong decision requires `|A_N-mu_j| >= c_*alpha_N/2`. NOC-003 gives

`P_piD,j(wrong) <= C[1/(M_N alpha_N^2) + N^-1/4/alpha_N + k_N^-1/2/alpha_N + (N/(sqrt(k_N)alpha_N))*pi_D,j(B_D)]`.

For the exact floored sequence, the three polynomial terms have orders

`1/(M_N alpha_N^2) = O(N^-1/8)`,
`N^-1/4/alpha_N=O(N^-3/16)`,
`k_N^-1/2/alpha_N=O(N^-13/32)`.

The fourth term decays faster than any power, by MT's positive departure-boundary energy and polynomial normalization. Hence the stationary reflected classifier error tends to zero. This is a **one-sided upper bound for a specific cosine test**, not an optimal risk rate or a practical N0.

Transfer to the actual unreflected iid-uniform labelled FIN by coupling the complete decision event on the observation interval. Under the inherited MT local mixing and exit bounds for each frozen operator, add

`delta_N,j <= mu0_N(D_N^c)+TV(nu_N P^D_{t_N},pi_D,j)+P_iid(τ_D<=t_N+M_N)`.

Each term is o(1) for a common sufficiently late polynomial deterministic origin `t_N` and polynomial `M_N`, both lying inside the accepted exponential residence window. Crucially the event-level TV error is **added once**, not divided by alpha_N and not multiplied by M_N. Thus

`P_actual,j(wrong) <= P_piD,j(wrong)+delta_N,j ->0`, for `j=0,1`.

**Disposition: CONDITIONAL_NATIVE_FINITE_N_SMOOTH_CLASSIFIER_CONSISTENCY.** The result inherits analytic MT/KO/DKO proof gates; the prior user-supplied DKO intake is an acceptance premise, not independently witnessed repository history. No finite-N lower/necessary resource bound is asserted; no shrinking-block sign-only theorem or universal operator identification follows.
