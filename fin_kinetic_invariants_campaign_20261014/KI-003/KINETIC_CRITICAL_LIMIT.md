# KI-003: prepared local kinetic-invariance theorem for fixed FIN A7

**Proof scope.** This is an analytic **conditional mathematical theorem** using the already independently accepted isolated fold (g_f,p_f,v), its one-null/ten-positive tangent Hessian, constant-diagonal A, and the existing microscopic count-state rules, not a new measured physical law. The source SF boundary completion is reused as a stopped martingale argument, not replaced by finite-dimensional convergence. All constants in estimates are existential, fixed independently of N, and valid after restricting to a sufficiently small interior fold neighbourhood.

## Statement
Fix lambda in a compact set, soft radius R>0, and a finite scaled time S. Let g_N=g_f+lambda N^(-2/3). For each K∈{HB,MET,BAR}, take the **exact** N-copy count generator with its frozen off-diagonal rates in KI-002. Prepare deterministic type states by balanced rounding of p_f+N^(-1/3)z0 r, with bounded |z0|<R, and no fitted kinetic auxiliary variables. For a sequence of such deterministic z0 or a stated tight mixture supported at distance >0 from the soft boundary with uniformly controlled initial fourth moments, define

x_K=ell_K^T(p_N-p_f), Z_{K,N}(s)=N^(1/3)x_K(N^(1/3)s).

Stop Z_{K,N} on reaching |Z|≥R or the fixed scaled horizon S. Then under the above hypotheses the stopped path converges in Skorokhod J1 to the unique stopped solution

 dZ_K=-m_K[a lambda+(b/2) Z_K^2] ds + sqrt(2m_K) dB_s.

For this **rounded initial class**, the *same* limit holds for the thermal/source-observation variable

 Z_{theta,N}(s)=N^(1/3)v^T[g_N X^T p_N(N^(1/3)s)-theta_f],

under the same stopping and absorption convention. For arbitrary initial states controlled only by E||w(0)||²<=C/N, the additional maximal transfer argument requires an appropriate fourth-moment/exceptional-event condition; it is not claimed for that larger class without it. Exit-side and bounded-horizon survival laws transfer as in SF-003's already admitted boundary completion; no uniform finite-N TV error or global exit law is asserted.

## Proof I: local limiting drift, including Metropolis kinks
On T use p=p_f+x r+w, x=ell_K^T(p-p_f), ell_K^T w=0. Since K_f is tangent-positive and H_F has precisely one null, B_K=(J_f restricted to ker(ell_K)∩T) is Hurwitz. Choose a fixed basis E_K on that invariant subspace, and P_K>0 solving B_K^T P_K+P_K B_K=-I. This P_K is only an analytical energy form, never a physical metric.

The finite-N count drift b_{N,K} differs uniformly from b_K(p,g)=-K(p,g)grad F_g(p) by O(N^-1). The field difference g A e_i/N is O(N^-1); all three attempt acceptances are globally Lipschitz on any compact gain/field range (MET's min is Lipschitz) and q is fixed.

At p_f, a tangent-space Taylor expansion of **smooth F**, and only local Lipschitz continuity of K, gives

 grad F_g(p) = const*1 + H_F(xr+w)+ (g-g_f)partial_g grad F_f +(x²/2)D³ F_f[r,r,·] + E,
 ||E||_tangent<=C(|x|³+|x|||w||+||w||²+|g-g_f|(|x|+||w||)),
 ||K(p,g)-K_f||<=C(|x|+||w||+|g-g_f|).

It is important that H_F r=constant*1 and K annihilates constants. Therefore the nonsmooth MET factor K(p)-K_f always multiplies a force of order |w|+x²+|g-g_f|; it cannot contribute a new O(x²) soft term. Applying ell_K^T K_f=m_K r^T and the identities a_F=a, b_F=b yields

 ell_K^T b_{N,K}(p,g_N) = -m_K[a(g_N-g_f)+(b/2)x²] + R_{N,K},
 |R_{N,K}|<=C(|x|||w||+||w||²+|x|³+|g_N-g_f|(|x|+||w||)+|g_N-g_f|²+N^-1).

Projection onto the fast invariant subspace similarly gives

 Q_K b_{N,K}(p,g_N)=B_K w+u_N,
 ||u_N||<=C(x²+|g_N-g_f|+(|x|+||w||)||w||+N^-1).

There is no assumed HB mean closure for MET/BAR.

## Proof II: stopped second-moment energy and no fast-tube exit
Let V=w^T P_K w (using coordinates of the chosen 10-dimensional basis). At a transition i→j,
 Δw=Q_K(e_j-e_i)/N, ||Δw||<=C/N, total off-diagonal rate <=C N.

The **exact finite-jump quadratic expansion** L_N V=2w^T P_K Q_K b_{N,K}+sum_{i,j} n_i r_ij(Δw)^T P_K Δw and Young's inequality yield for a sufficiently small fixed interior neighbourhood

 L_N V <= -k V+C(x⁴+|g_N-g_f|²+1/N).

For stopping at |Z|≥R or ||w||≥N^-1/4 or T_N=N^(1/3) S, the nonnegative stopped terminal V is kept in Dynkin's formula. With balanced rounded initial w(0)=O(N^-1) it follows that

 E∫_0^(T_N∧stop) ||w_t||² dt=O(N^-2/3),
 P(fast tube first) =O(N^-1/6).

As in SF, the expected integrated rescaled soft drift errors have negative exponents: x w gives N^-1/6; w² N^-1/3; x³ and delta_g*x give N^-1/3; delta_g*w gives N^-1/2; delta_g² gives N^-2/3; N^-1 gives N^-1/3. The noise covariance at the stopped occupations converges to D_f because ||p-p_f||<=C_R N^-1/3+N^-1/4+O(N^-1). No claim about stationary preparation is needed.

## Proof III: a fourth-moment maximal estimate (required for same-observable transfer)
For the **same** energy V, the exact generator identity for its square is

 L_N(V²)=2 V L_N V +sum_{i,j} n_i r_ij [V(w+Δw)-V(w)]².

Since |ΔV|<=C(sqrt(V)/N+N^-2), the positive jump contribution is <=C(V/N+N^-3). Combining with the previously paid LV inequality and Young's inequality gives a fixed local inequality

 L_N(V²)<=-k' V²+C(x^8+|g_N-g_f|^4+N^-2).

This uses **actual bounded jumps** and a uniform rate O(N), not an assertion that the stable modes are already at stationarity. While |x|<=C_R N^-1/3, the forcing is O(N^-2). Stopped Dynkin, with the nonnegative V² at the terminal exit, and V(0)²=O(N^-4), implies

 E V²(w_{T_N∧stop}) <= C N^-5/3.

For any fixed epsilon>0, stop instead at the *earlier* first crossing of ||w||≥epsilon N^-1/3. This threshold is eventually smaller than N^-1/4 and lies inside the local proof domain. Since V²>=c epsilon^4 N^-4/3 upon crossing (including overshoot),

 P(sup_{t<=T_N∧soft-exit} ||w_t|| >= epsilon N^-1/3) <= C epsilon^-4 N^-1/3 + P(original auxiliary exit).

In fact the earlier-stop argument absorbs the auxiliary event; in either form the right-hand side tends to zero. Thus, on compact soft-stopped windows,

 sup ||w_t||=o_P(N^-1/3).

This stronger statement follows for the frozen rounded preparation; the previous SF second-moment estimate alone would **not** have yielded this comparison.

## Proof IV: common thermal readout
Both ell_HB=g_f Xv (modulo constants) and ell_K satisfy ell*r=1. Hence

 ell_HB^T(p-p_f) = x_K + (ell_HB-ell_K)^T w.

Also the exact finite-N coordinate conversion is

 v^T[g_N X^T p-theta_f] = (g_N/g_f)ell_HB^T(p-p_f) + [(g_N-g_f)/g_f] v^T theta_f.

Multiplying by N^(1/3), the deterministic gain difference is O(N^-1/3), and the stochastic discrepancy is N^(1/3)O_P(sup||w||)=o_P(1). The two readouts therefore have the **same compact-stopped weak limit** under the rounded initial class, even though their finite-N observables differ and MET/BAR do not have exactly closed mean dynamics.

## Proof V: compensated martingale and absorbing boundary
Each jump of Z_{K,N} has size O(N^-2/3). On accelerated time the compensator's integrated drift converges in L1 to -m_K[a lambda+b Z²/2] by Proof II, and its predictable bracket converges uniformly before stopping to 2m_K ds because D_f=2K_f and full occupation approaches p_f. The integrated absolute third-jump characteristic is O(N^-2/3). Bounded drift, tight bracket and vanishing jump sizes imply tightness on stopped D([0,S]) with continuous limits. Any limit solves the bounded-domain martingale problem for G_K phi=-m_K(a lambda+b z²/2)phi'+m_K phi''. Smooth drift and strictly positive constant diffusion give stopped uniqueness.

For **boundary absorption**, reuse SF-003's accepted joint convergence method: track (stopped coordinate, stopped clock, stopped martingale) in a common joint filtration; work with an outer radius R+epsilon; nondegenerate Brownian crossing makes the crossing map continuous at almost every limiting path. Capped exit times and upper/lower/survival records follow. The vanishing auxiliary-tube exit is handled before removing the auxiliary stopping condition. Merely checking interior test functions would NOT establish this boundary claim.

## Conventions and exclusion
The theorem is dimensionless, prepared, local and compact-stopped. It gives no finite-N rate with certified constants, no global mixing or B4, no physical seconds, no irreversible physical bath, and no claim for arbitrary reversible kinetic rules or arbitrary initial laws. The exact fold root and inertia enter as already accepted interval-assisted premises; all specific decimal mobility values are binary64 proposals.
