# RI-001 — Exact tagged/count parent and disjoint observation contracts

**Scientific disposition:** exact algebra for the supplied labelled LOO-HB parent; microscopic transfer uses admitted OBS/MT intake as a *conditional dependency*. No extra transition law.

Let `xi=(xi_a)_{a=1}^N`, `n_i=#{a:xi_a=i}`, `b=(xi_{a_1},...,xi_{a_k})`, `c_i(b)=#{h:b_h=i}`, and `R_ij(n)=softmax_j[(g_N/N)A(n-e_i)]`. Valid states are `n_i>=c_i(b)` and `sum_i n_i=N`. For i≠j the `(n,b)` generator has:

1. Tagged transition h with b_h=i: `(n,b)->(n-e_i+e_j,b^(h->j))` at rate `R_ij(n)` per tagged coordinate.
2. Untagged transition i->j: `(n,b)->(n-e_i+e_j,b)` at rate `(n_i-c_i(b))*R_ij(n)`.
3. Diagonal equals the negative sum of all genuine off-diagonal rates. Self-refresh `j=i` changes no state, has zero generator contribution, but is not a new attempt clock. The per-labelled-copy attempt clock remains 1.

This is an exact Markov quotient of the supplied labelled generator: on functions of `(n,b)` all same-fibre labelled states have the same transitions; summing over tagged and untagged departures recovers `n_i R_ij(n)`. The **count** process alone does not determine tagged histories unless its labelled generator is also declared.

Distinct sigma fields, for each fixed critical-time interval `[0,S]` and tags frozen independently of the path:
- `TAG_PARITY_PATH_N`: `sigma(S_{a_h}(t_N+sqrt(N)s), h<=k, s<=S)`.
- `TAG_Q12_PATH_N`: same with the exact q12 label xi rather than parity; strictly richer in general.
- `EXACT_INTEGRAL_PATH_N`: `sigma(Y_N(s),s<=S)` with `Y_N=N^(-1/4) int_{t_N}^{t_N+sqrt(N)s}m_B(t)dt`.
- `QUANTIZED_TWO_INCREMENT_N`: sigma field of the pair `(Q(V1),Q(V2))` for the fixed RI-002 settings.
- `CONTINUOUS_MESOSCOPIC_Y`: `sigma(Y_s,s<=S)` in the limiting diffusion; not assumed to have the same information as the exact microscopic integral.

Since `dY_N/ds=N^(1/4)m_B(t_N+sqrt(N)s)` wherever the derivative exists, right derivatives recover every piecewise-constant block average except its jump times. **For k=1** the whole exact integral path recovers the parity history, modulo null time sets and choice of right-continuous versions. This inverse is not continuous under weak path convergence. For k>1, recovering only the *block mean* is weaker than recovering each tagged path.

Example of lost inference under weak convergence: C is fair ±1; record `T_N=C/N`; then `T_N->0` uniformly yet sign(T_N) recovers C perfectly for every N. This is **not** a FIN counterexample; it refutes the generic inference from weak path convergence to full-record sufficiency.

Small exact support checks in `independent_checks.py` include the required `(N,M,k1,k2,h1,h2)=(12,2,4,3,3,0)`, whose joint probability is exactly zero.

The original labelled model and MT iid-start window remain external premisses (source hashes in `SOURCE_HASHES.json`); this package does not include or modify the previous archives.
