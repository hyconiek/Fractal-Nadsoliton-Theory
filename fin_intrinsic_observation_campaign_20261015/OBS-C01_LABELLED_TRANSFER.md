# OBS-C01 — Exact labelled parent, orbit lift and local-window transfer

**Scientific status:** exact labelled/count identities; scoped local-stationary path transfer conditional on accepted MT/UF analytic premises. No new microscopic update rule.

## Frozen identities

Let q=12, sigma_i=(-1)^i, n_i=# {a: xi_a=i}, N in 12ℕ and g_N=g_c+lambda N^{-1/2} with **exact implicit** g_c=12/Lambda6. Unit-attempt leave-one-out heat bath is

L^lab f(xi)=sum_{a=1}^N sum_{j=0}^{11} s_j((g_N/N) A(n-e_{xi_a}))[f(xi^{a->j})-f(xi)].

No per-slot permutation moves occur. The map xi->n intertwines L^lab with the exact MT count generator. For every permutation pi of the N slots, L^lab(f∘pi)=(L^lab f)∘pi. The iid uniform initial measure is slot exchangeable; covariance of the semigroup preserves exchangeability at each deterministic time. The finite orbit of configurations having the same exact n consists of N!/∏n_i! elements. Conditional on n the labelled law is uniform on that orbit. This is true at every deterministic time, NOT an independence assertion across times.

Define the Markov **orbit lift** H(ν)=sum_n ν(n) Unif{xi: counts(xi)=n}. Count fibres are disjoint, so TV(Hν,Hη)=TV(ν,η); the projection and any fixed tagged observation contract TV. Therefore accepted MT count-window bound epsilon_N(t)=TV(Law_iid(n(t)),pi_D,N)->0 uniformly for t_loc(N)<=t<=exp(cN) transfers exactly, *at a single time*, to the full labelled orbit lift H(pi_D,N).

## Path transfer from the selected time

Freeze deterministic t_N satisfying t_N>=t_loc(N) and t_N+sqrt(N)S<=exp(cN) for each desired fixed S and fixed c<b. At t_N, couple the entire labelled vector with H(pi_D,N); the chance of discrepancy is epsilon_N(t_N). Subsequently run both **the same labelled HB kernel**: total variation of their complete future-path distributions is no larger than epsilon_N(t_N). H(pi_D,N) is reversible/stationary for the labelled chain with *only count-outgoing jumps suppressed*, because the finite parent is detailed balanced against the all-pairs labelled Gibbs measure and the conditional count weights include multinomial factors.

The reflection is analytical: its internal edges are unchanged and it can be coupled exactly with the unreflected parent up to the first D_N exit. From H(pi_D,N), the probability of a suppressed boundary transition before sqrt(N)S is bounded by sqrt(N)S times a polynomial prefactor times exp(-Nb), using MT-001 boundary height and local partition lower bound. It tends to zero. Use this additional term; a count-path MT window alone does not let one erase exits.

The conditional MT local Gibbs has the Z=N^(1/4)m marginal converging to rho_lambda(z)∝exp(chi lambda z²/2−z⁴/12), and E||w||²=O(1/N), with higher polynomial moments and tails, where w=p−p_m. The UF stopped generator proof extends from iid-at-time0 to this random initial law because its hypotheses need (i) Z(0) tight and convergent, (ii) E||w(0)||²=O(1/N), (iii) auxiliary tube exceptional probability vanishing, (iv) full interior support asymptotically; all are paid by MT-003. The small initial critical Z law need not be a delta at 0. Localisation plus nonexplosion extends convergence on any fixed scaled horizon to the quartic diffusion initialized at **rho_lambda**, making the limit a stationary continuous process (stationarity follows from the 1D generator, not from exchanging N and infinite-time limits).

This establishes comparison in distribution of the joint *count path*, and, by the frozen labelled generator, licences subsequent tagged-history martingale analysis. It does **not** license tagging arbitrary count processes: slot histories are extra model content.

## Outside-primary negative control

Artificially swapping labels between two slots at rate1 changes labelled-tag histories while preserving every count path. Such a generator is labelled OUTSIDE_PRIMARY_FIN; it is not a counterexample to the fixed labelled HB law. Tags are labels, not a graph metric or spatially adjacent physical subsystems. No source of physical data/time/hardware is inferred.
