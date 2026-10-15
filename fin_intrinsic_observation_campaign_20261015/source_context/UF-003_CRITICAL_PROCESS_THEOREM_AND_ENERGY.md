# UF-003 — Prepared-reference-start critical process theorem

## Theorem — iid-uniform critical path limit (conditional on admitted exact input spectrum)

Assume the exact-decimal q12 A with the **accepted unique simple top k6 eigenvalue** and a strictly negative ten-mode restriction at g_c=12/lambda6. Fix lambda∈R, S<infinity. For each N start N independently iid uniform labels, evolve the actual rate-one-per-copy leave-one-out HB generator at g_N=g_c+lambda N^-1/2, and define Z_N(s)=N^1/4 sigma^T p_N(N^1/2 s). Then in Skorokhod J1 on [0,S],

 Z_N => Z solving dZ=(chi lambda Z-Z^3/3)ds+sqrt(2)dB, Z(0)=0.

The unique solution is nonexplosive. The same result holds separately for the two **explicitly declared** g0-matched MET/BAR diagnostic generators. It is NOT a statement about general reversible controllers or full Gibbs initial laws. For every fixed R>0, capped first hitting times and exit side/survival outcome for (-R,R) converge as well. The theorem uses earlier input spectral acceptance and ordinary analytical smoothness/Lipschitz facts; it is not mechanized or independently peer-reviewed.

## Proof 1: ten-mode local stability and actual generator inequality

Work with orthogonal fast projector P_f=I-11^T/12-sigma sigma^T/12 and p=p_m+w. Exactly b_fast(p_m,g)=0 for each rule by shift-two covariance. At (u,g_c), all three rules share K(u,g_c)=Sigma_u and derivative J=-Sigma_u Hess F_{g_c}, which equals (g_c/12)A-I on tangent. On the ten-dimensional fast space its eigenvalues are strictly negative by the supplied accepted spectral input. This J is **symmetric** in the ordinary Euclidean fast metric, so the analysis energy V=||w||^2 needs no invented metric, singular inverse or finite-N spectral-gap claim.

Because K is locally Lipschitz and grad F is smooth on a fixed interior neighbourhood, the fast drift satisfies b_fast(p_m+w,g)=J_f w+O((|m|+|g-g_c|)||w||+||w||^2). The actual finite-N LOO drift deviates by O(1/N). Each jump changes w by P_f(e_j-e_i)/N with norm <=sqrt(2)/N, and total rate is <=C N. For V=||w||^2 the exact quadratic generator increment is

 L_N V = 2w^T P_f b_N(p,g)+ sum_{i!=j} n_i R_ij^(N) ||P_f(e_j-e_i)/N||^2.

The last term is nonnegative and <=C/N (never discarded). Shrink a **fixed** interior neighbourhood, not an N-dependent stability assumption, and apply Young's inequality. For large enough N and every state stopped in this neighbourhood,

 L_N V <= -k V + C/N,

for positive constants k,C independent of N. No claimed numerical k,N0 is needed for an asymptotic theorem. The same argument is valid for MET's locally Lipschitz K; the differentiability needed at u follows because gradF(u)=const·1 and K1=0, not because Metropolis rates are C3.

## Proof 2: unbiased initial moments, energy occupation and fast tube

Take soft stop |Z_N|>=R, fast stop ||w||>=delta_N=N^-1/5, and physical model horizon T_N=N^1/2 S. The initial iid moments E V(0)=O(1/N), E V(0)^2=O(1/N^2), and initial bad-event probability o(1) were paid in UF-001. Stopped Dynkin's formula, retaining the **nonnegative terminal energy** and the positive jump overshoot, gives

 E V(T_N∧tau)+k E ∫_0^(T_N∧tau) ||w||² dt <= E V(0)+C T_N/N =O(N^-1/2).

Consequently

 E ∫_0^(S∧tau/N^1/2) ||w(N^1/2s)||² ds = O(1/N).

If the fast tube is reached before the soft stop/horizon, the nonnegative terminal V is at least delta_N² (the jump overshoot only increases it). Therefore

 P(fast-tube-first) <= C N^-1/2 /delta_N²=O(N^-1/10).

This includes the initial tail event by the same inequality. One must not write E V(t∧tau)<=C e^-kt+C/N: after absorption the exit value is not subject to the original contraction. Instead the above stopped integral and terminal inequalities are used.

## Proof 3: pay each dangerous term and identify the drift

On the soft stop, |m|<=R N^-1/4+2/N (jump overshoot), |delta g|<=C N^-1/2. Use the exact selection-rule estimate from UF-002. In scaled time the drift-multiplier is N^(3/4). The following integrated terms vanish in L1:

* N^(3/4) E ∫||w||² ds = O(N^-1/4).
* N^(3/4)(|delta g||m|+|m|³) E∫||w|| ds <= C N^(3/4) N^-3/4 O(N^-1/2)=O(N^-1/2), by Cauchy–Schwarz.
* N^(3/4)|m|^4=O(N^-1/4).
* N^(3/4)|delta g|m²=O(N^-1/4).
* N^(3/4)|delta g|²|m|=O(N^-1/2).
* N^(3/4)*O(1/N)=O(N^-1/4), the **actual LOO rate** payment.

The gain-linear term and cubic term converge exactly to chi lambda Z-Z³/3. If one uses stronger HB expansions some bounds improve but none is needed.

## Proof 4: predictable martingale bracket, jumps, tightness

The actual Z_N compensated martingale has accelerated predictable bracket density equal to B(p_N,g_N), where B=sum p_i R_ij^(N)(sigma_j-sigma_i)^2, with the correct clock/scaling cancellation N^(1/2)*N^(1/2)/N=1. Since |m|<=R N^-1/4 and ||w||<=N^-1/5 on the stop, B=2+o(1) uniformly; the bracket tends to 2(s∧tau_R). A jump of Z has size at most 2N^-3/4, and its accelerated absolute third characteristic is O(N^-3/4) on bounded horizons. Every compactly stopped drift integral has bounded variation in probability: the polynomial limit is bounded on [-R,R] and the integrated remainder is o_L1(1). These facts yield tightness (Aldous-Rebolledo) and identify limit points by the joint martingale problem. UI of bounded-on-stop quadratic martingale/bracket and the square-martingale identity follows by localization and bounded predictable brackets; the limiting martingale is a Brownian time-change in its joint limiting natural filtration, *not* by pretending the random stop is deterministic.

The killed/interior problem for the one-dimensional generator

 Gf(z)=(chi lambda z-z³/3) f'(z)+f''(z)

has a unique stopped law on (-R,R), because drift is locally Lipschitz and diffusion sqrt(2) is nonzero. To transfer the soft boundary, first work at larger radius R'>R, identify the unfrozen continuous process there, and invoke a.s. continuity of first hitting/crossing at ±R for nondegenerate one-dimensional diffusion; then handle the stopping/cemetery convention jointly with the martingale. The fast-stop probability is O(N^-1/10), so it can be removed in probability. Finite-horizon side/SURVIVE outcomes are continuity sets (no atom at a fixed S).

## Proof 5: remove the soft stop for finite horizons

The limiting diffusion is nonexplosive. Indeed G(z²)=2chi lambda z²-(2/3)z⁴+2<=C_lambda for all z; optional stopping gives P(sup_{s<=S}|Z_s|>=R)<= (E Z0²+C_lambda S)/R². Apply compact-stopped convergence at every fixed R and then R->infinity. This diagonal localization establishes full D([0,S]) weak convergence. It does **not** give convergence uniformly as S->infinity or the finite-N equilibrium process.

## Scope and dependence

This is an analytical proof submitted for independent mathematical intake, conditional on correct previously admitted exact spectral binding. It invokes standard stopped semimartingale tightness/uniqueness machinery with the bounded-bracket/UI argument detailed above; it is not formal proof-assistant verification or experimental validation. No explicit numeric finite-N TV bounds, physical clock, source of kernel, all-preparation theorem or global metastable exit is concluded. The MET control's kink is paid via Lipschitz mobility and never replaced by a fake C3 derivative.
