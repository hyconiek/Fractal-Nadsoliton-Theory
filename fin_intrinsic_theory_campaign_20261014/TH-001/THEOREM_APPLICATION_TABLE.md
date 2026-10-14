# TH-001 theorem application table

## L01 — Kraaij (2015), Theorem 2.3: directly applied in TH-002

For the q-state Glauber-type empirical measure chain, Theorem 2.3 assumes a continuous limiting edge intensity `v(a,b,p)`, uniform convergence of `p_a r_N(a,b,p)` to that intensity, positivity whenever `p_a>0`, a monotone boundary factorization when `p_a=0`, and an initial-law LDP.

For FIN define

`v_ij(p)=p_i s_j(p)`, `s(p)=softmax(g A p)`.

Condition-by-condition:

1. **Continuity:** PASS. `s_j` is analytic on the compact simplex.
2. **Interior positivity:** PASS. `s_j(p)>0` for all finite logits, hence `v_ij>0` whenever `p_i>0`.
3. **Boundary factorization:** PASS with `v_dagger(x)=x` and `v_ddagger(p)=s_j(p)`. The first factor is increasing; the second is continuous and strictly positive at the boundary.
4. **Uniform rate convergence:** PASS for each fixed finite `g`, and uniformly on compact gain intervals. The finite rate is `p_i softmax_j[gA(p-e_i/N)]`. Coordinate softmax has l-infinity Lipschitz constant at most 1/2, so the discrepancy is at most `g max|A_ki|/(2N)`.
5. **Initial LDP:** CONDITIONAL in the theorem, but two concrete classes are paid here: deterministic type sequences converging to `p0`; and the declared Gibbs equilibrium sequence, whose type probabilities obey the method-of-types LDP with rate `F_g-min F_g`.

Therefore TH-002 may invoke Theorem 2.3 for these initial classes. This is an application to an interacting mean-field chain, not the independent-walker result ST8548.

## L02 — Mielke–Peletier–Renger (2014)

Used as structural context only. Their paper connects path large-deviation costs with generalized gradient structures. TH-002 does **not** transfer an independent-particle theorem to interacting FIN. Instead, the FIN detailed-balance/reversal identity is derived directly from its Hamiltonian and free energy.

## L03 — Cuff et al. (2012)

Benchmark only. The paper proves an `n^(-2/3)` spinodal window and `n^(4/3)` mixing scale in single-update steps for mean-field Potts. Since FIN uses total attempt rate N per model-time unit, `n^(4/3)` single updates corresponds dimensionally to `n^(1/3)` FIN model time. No Potts coefficient or mixing theorem is transferred.

## L10 — Katzenberger (1991)

Not invoked as a completed theorem application. It is a possible method for strongly stable directions, but its SDE/manifold hypotheses are not automatically a theorem for the present density-dependent jump chain. TH-005 isolates the exact fast-mode/perturbed-test-function lemma that would have to be paid before using such a method.

## Scope boundary L04

Uniform-in-time mean-field convergence results requiring one exponentially stable nonlinear stationary law are not used at the fold: the fold has a zero mode and nearby multiple stationary branches, so that premise is unavailable exactly where TH-005 operates.
