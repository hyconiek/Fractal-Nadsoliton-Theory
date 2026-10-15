# DKO-001 — Characteristics derived from the labelled FIN generator

Date: 2026-10-15. **Analytic scope:** fixed fast model-time horizon, q12 unit-attempt leave-one-out heat bath, N in 12N, g=g_c, two frozen model IDs ε=0,+1/100. These are microscopic generator calculations, not assumed OU rates.

## 1. Exact, finite-N formulas

A labelled configuration x=(x_a)_{a=1}^N has population n_j=number of labels j and a deterministic block B of k slots. Set α_N=k/N, c_a=N/k−1 for a in B and c_a=−1 otherwise, f_r(j)=√2 cos(2πrj/12), r=1,3,4. Put

F_{r,N}=N^{-1/2}Σ_a f_r(x_a),  U_{r,N}=N^{-1/2}Σ_a c_a f_r(x_a), V=F+U.

For a slot a with current label j, its refresh to l has rate q_l(n−e_j)=softmax[(g/N)A(n−e_j)]_l. Null refreshes have zero state increment, and must not be counted as genuine jumps. Let Δ_{jl} f=f(l)−f(j), and a_f(v)=Σ_l f(l) softmax[(g/N) A v]_l. Directly from the ORIGINAL generator:

L F_r=N^{-1/2}Σ_a [a_{f_r}(n−e_{x_a})−f_r(x_a)],
L U_r=−U_r+N^{-1/2}Σ_a c_a[a_{f_r}(n−e_{x_a})−a_{f_r}(n)].

The latter is **exact** (Σ_a c_a=0); each term in the remainder is O(1/N) uniformly in configuration for bounded g,A,f. Thus |L U_r+U_r|≤C_α/√N for α_N bounded away from zero. This pays the private decay; it is not assumed from static covariance.

For X_r∈{F_r,U_r}, the exact predictable bracket density is, with weights d_a(F)=1,d_a(U)=c_a,

d⟨M_{X_r},M_{X_s}⟩/dt=(1/N)Σ_a d_a(X)d_a(Y) Σ_l q_l(n−e_{x_a}) Δ_{x_a l}f_r Δ_{x_a l}f_s.

For r=s, the inner conditional square is E_q[f_r(new)^2]+f_r(old)^2−2f_r(old)E_q[f_r(new)]. At uniform it averages to **2**, but for individual old labels it is 1+f_r(old)^2, NOT universally 2. In particular the parity-only shortcut fails for a general cosine harmonic.

## 2. Fast-global drift and dangerous soft terms

Let p=n/N=u+(m/12)σ+w, where σ_j=(−1)^j, σ·w=1·w=0. The softmax Lipschitz bound gives

LF_r=√N[f_r·softmax(gAp)−f_r·p]+O(N^{-1/2}).

On the exact parity line p_m, f_r·softmax(gAp_m)=f_r·p_m=0 because r=1,3,4 are transverse to {1,σ}. At u, the derivative along mode f_r is (gΛ_r/12−1)=−γ_r, with γ_1=1 and γ_{3,4}=1−Λ_{3,4}/Λ_6. Smoothness on the inherited small positive cylinder and parity covariance yield

|LF_r+γ_r F_r|≤C[|m|√N||w||+√N||w||²+N^{-1/2}]

(with constants uniform for the two frozen operators). No order-one pure-soft forcing survives. The parity coordinate need not be artificially held at zero; m=O_P(N^{-1/4}). These are valid only with local occupation moments from MT; a bare static snapshot theorem is insufficient.

## 3. Brackets and block sampling

The exact bracket above, expanded at the uniform population and using the conditional block hypergeometric law, gives, in integrated probability on every fixed fast horizon T,

⟨M_{F_r},M_{F_s}⟩_T→2δ_rs T,
⟨M_{U_r},M_{U_s}⟩_T→2c_α δ_rs T,
⟨M_{F_r},M_{U_s}⟩_T→0,

where c_α=(1−α)/α. Indeed (1/N)Σ_a c_a=0 and (1/N)Σ_a c_a²→c_α. Also the uniform covariance E_u f_r f_s=δ_rs. Conditional uniformity within count orbits gives E||p_B−p||=O(N^{-1/2}) for fixed α; MT moments give E||p−u||=O(N^{-1/4}). Hence the mean absolute error of every bracket integrand vanishes. The same argument handles different Fourier coordinates. Here *integrated probability* suffices; no false supremum-in-time high-probability fast norm is asserted.

## 4. Initialization, transfer, and bounded jumps

MT's conditional Gibbs joint quartic/fast Gaussian limit and KO's **static** conditional hypergeometric CLT give, jointly at the common origin: soft Z distributed as quartic ρ, F_r∼N(0,γ_r^{-1}) and U_r∼N(0,c_α), mutually independent across the listed coordinates and independent of Z in the limit. Initial independence is asymptotic only. Every F/U jump is bounded by C_α/√N, with total refresh attempt rate N, hence the limiting paths are continuous.

The stationary reflected *labelled* lift is used solely in the proof. Its one-time expectations satisfy MT moment bounds uniformly in time; integrated moment/bracket errors on [0,T] consequently vanish by Markov's inequality. A path coupling transfers the result to the actual iid-start parent at deterministic post-MT t_N, retaining initial outside-D and exit probabilities. It does NOT move unbounded moments through o(1) path TV.

**Scope:** all conclusions conditional on the accepted MT/KO local domain, mixing and labelled transfer for both frozen operators. No original review bytes are invented from the user-pasted intake. No physical time unit, spatial Fourier wave or new microscopic law.
