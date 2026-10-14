# UF-004 — Local critical law versus global finite-N Gibbs

## The derived limiting diffusion invariant law

For the accepted UF-003 finite-horizon path theorem, the one-dimensional generator is G f=(chi lambda z-z³/3)f'+f''. Solving the zero-current equation rho'/rho=(chi lambda z-z³/3) gives the **normalized diffusion invariant density**

 rho_lambda(z)=Z(lambda)^(-1) exp[chi lambda z²/2-z⁴/12],
 Z(lambda)=∫_{-∞}^∞ exp[chi lambda z²/2-z⁴/12]dz <∞.

A Lyapunov condition and nondegenerate 1D noise give uniqueness/ergodicity for this mathematical limiting diffusion. At lambda=0 its even moments have the exact Gamma expression

 E|Z|^k = 12^(k/4) Gamma((k+1)/4)/Gamma(1/4), k>-1.

In particular E Z²=sqrt(12) Gamma(3/4)/Gamma(1/4), nominally about 1.17082865661. This is not a full finite-N Gibbs equilibrium prediction and not a fundamental FIN constant. The limiting density is symmetric: there is **no strict sign-selection rule**. A symmetric microscopic initial distribution is invariant under cyclic label shift by one, which sends m to -m; the trajectory readout law is symmetric at every N and for each of the three shift-equivariant kinetics.

## Direct independent exact-model free-energy competitor

With F_g(p)=sum_j p_j log(12p_j)-(g/2)p^T A p, F_gc(u)=0 because Au=0 and entropy relative to u vanishes. For a one-label occupation e_0, F_gc(e_0)=log12-(g_c/2)A_00.

Using only the admitted A7 construction and exact g_c (not a numerical fitted gain), the binary64 proposal gives F_gc(e_0)=-0.773026708896337. It is strictly negative with a large margin: using the earlier accepted spectral boxes (or coarse interval bounds) can prove that gc*A00/2>log12; the present package supplies no new outward interval certificate for that bound, so this particular arithmetic cross-check is nominal, supported independently by the earlier admitted coexistence/equal-energy competitor in the supplied intake.

More importantly, the **inherited admitted coexistence statement** at g_eq<g_c, with F_{g_eq}(p_eq)=0 and p_eq^T A p_eq>0, gives by exact algebra:

 F_{g_c}(p_eq)=-(g_c-g_eq)(p_eq^T A p_eq)/2<0.

This establishes a strict global free-energy competitor independent of nominal e_0 rounding. No first-barrier Gamma/capacity theorem is assumed. Since F_gc is continuous and F_gc(u)=0, a fixed sufficiently small neighbourhood U of u has inf_{p∈U} F_gc(p) greater than the competing negative value. The finite-N Gibbs method-of-types upper/lower bounds imply P_Gibbs(p_N ∈ U)→0 exponentially as N grows for gains g_N→g_c. An interior competitor close to p_eq can be used if desired. This is a theorem about the **full occupation vector**, not a licensed assertion about the global equilibrium distribution of the scalar parity m alone.

## Distinct limit orders

1. For iid-uniform microscopic initial labels, first N→∞ at each fixed scaled time s (time=N^1/2 s), then s→∞: the full occupation converges in probability to u for every fixed s (m=N^-1/4 Z_N and w→0 on fixed horizons), while the magnified parity variable Z_N converges to the quartic diffusion. That limiting diffusion subsequently relaxes to rho_lambda as s→∞.

2. For each fixed N, first physical model time t→∞, using the finite irreducible reversible Gibbs chain; then N→∞ with g_N→g_c: full occupation lies outside some fixed neighbourhood of u with probability tending to one. This is inconsistent with the full-occupation u limit in order 1.

The upper/lower law limits do not commute **for the explicitly named full-occupation observable**. No automatic scalar-parity difference is asserted, since ordered equilibrium branches may project to the same parity value. No quantitative global switching time, no uniform-in-time approximation and no limiting full-parent stationary quartic law follow from UF-003.

## Scientific status

Diffusion invariant density and its moments: **exact mathematical consequence of proved limiting SDE**. Gibbs non-dominance near u: **conditional on the inherited exact coexistence claim; independently corroborated by a strong numeric one-label competitor**. The local dynamical process does not describe the entire Gibbs ensemble. Neither unbiased mathematical initialization nor tuning g_N to g_c is a physical origin story or unique sign selector.
