# MT-001 — Local coercivity, one connected critical type domain

**Scientific status:** analytic proof conditional on the already accepted **exact** FIN spectral data: q=12, A=XX^T positive semidefinite/circulant, top parity eigenvalue Λ6 simple, g_c=12/Λ6 in (5,6), and fast relative gap γ=1−max(Λ3,Λ4,Λ5)/Λ6 > 0.0186. These are inherited model premises, not independently reconstructed rigorous intervals here. N ranges over positive multiples of 12; |λ|≤Λ finite. All constants N0 can depend on Λ, but *not* on N.

## A. Explicit analysis region

Let u=1/12, σ_j=(-1)^j, m=σ·p, p_m=u+mσ/12, w=p−p_m (w⊥1,σ). Set ε=1/100, δ=ε²=1/10000, h=1/10^10, b=h/2 and

 C={p: |m|<ε, ||w||_2<δ, p_i>0},  D_N=the nearest-neighbour count-graph component containing n_*=(N/12,...,N/12) of {n/N∈C, F_{g_N}(n/N)≤h},
 g_N=g_c+λ N^(-1/2), F_g(p)=Σ_i p_i log(12p_i)−g p^T A p/2.

All adjacency uses **actual single-copy count edges** n→n−e_i+e_j (i≠j); the generator/rates are not changed. Strict cylinder inequalities and the closed sublevel are exact domain definitions. The outer boundary is artificial only for analysis.

## B. Exact parity and transverse structure

Circulant A and parity imply A p_m=χ m σ, χ=Λ6/12=1/g_c. Since F's gradient at p_m is constant inside each parity class, its projection onto W={1,σ}⊥ vanishes: ∂_w F_g(p_m)=0 for every |m|<1, g. The exact one-dimensional free energy is

 F_g(p_m)=I(m)−(gχ/2)m²,  I(m)=((1+m)/2)log(1+m)+((1−m)/2)log(1−m).

For |m|<1, I(m)=m²/2+m⁴/12+m⁶/30+... (positive coefficients). Thus

 F_{g_N}(p_m)≥m⁴/12−χλ m²/(2√N)≥m⁴/24−(3/2)χ²Λ²/N.  (A)

Here the second bound completes the square in m².

For p in C, p_i ≤ (1+ε)/12+δ. The fast Hessian is the restriction of diag(1/p_i)−g_N A. The inherited gap, g_c>5, |g_N−g_c|≤.001 for N large and the rational estimates Λ_fast≤12/g_c<12/5 yield

 H_W≥[1/((1+ε)/12+δ)−12(1−0.0186)−0.001*(12/5)] I_W > 0.08 I_W.  (B)

Every numerical constant in (B) is rational except the inherited **strict spectral premise**, and can be checked with Fraction. Integrating this bound from p_m to p_m+w, with vanishing fast gradient at p_m, gives

 F_{g_N}(p_m+w) ≥ m⁴/24+0.04||w||²−C_Λ/N, C_Λ=3χ²Λ²/2. (C)

For sufficiently large N, (C) gives at |m|=ε a lower bound approaching 4.1666...×10^-10 and at ||w||=δ a lower bound approaching 4×10^-10. Both strictly exceed h=10^-10. The entire sublevel is therefore separated by a fixed positive distance from the outer cylinder. In particular the type steps (size √2/N) cannot leave C while remaining in the sublevel, once N≥N0.

## C. No hidden disconnected count component

Introduce a **central logarithmic box**

 C_N={|m|≤a_N, ||w||≤b_N}, a_N=N^(-1/4)(log N)^(1/4), b_N=N^(-1/2)(log N)^(1/2).

Both thresholds tend to zero and exceed the critical equilibrium scales by an increasing logarithmic factor. Its finite-label type points can be linked to n_* as follows. First redistribute counts *within each of the six even labels* and *within each of the six odd labels* until each parity class is balanced to within one count per label; this holds m fixed and keeps the deviation of each count from that parity mean bounded by its starting deviation plus O(1). Next transfer one count at a time between opposite parities toward m=0, rebalancing the two six-label classes after each transfer. Every intermediate count is nonnegative and is a legal single-copy edge. Along the constructed route, |m|≤a_N+O(1/N), ||w||≤C(b_N+1/N), and the length is ≤C_q N (fixed q). The local Taylor upper bound, uniform for |λ|≤Λ,

 |F_{g_N}(p)|≤C(m⁴+||w||²+Λ²/N)

holds on this route, so the energy stays ≤C log N/N<h for all sufficiently large N. Thus all central types lie in D_N. **No exact finite-N parity lumpability is used.**

To connect an arbitrary other type p in {C,F≤h}, use the actual gradient ∇_T F of the smooth *continuous* free energy only as a path-selection device. The strict fast Hessian and parity expansion imply on a sufficiently small local region (the chosen C is small enough using bounded third derivatives)

 ||∇_T F_{g_N}(p)||≥c(|m|³+||w||)−C|λ|N^(-1/2)|m|. (D)

The gradient payment can be made explicit at ε=.01, rather than relying on an unverified small-domain assertion. Let p_min=(1−ε)/12−δ=0.0824 and C_m=1/(24 p_min²)<6.14. Because the fast gradient vanishes at p_m and the fast Hessian is ≥.08 on the segment, ||∇_w F||≥.08||w||. Differentiating the integral Taylor remainder twice in w and once in m uses only the entropy Hessian (the −gA term has zero m derivative), giving |∂_m F(p_m+w)−∂_m F(p_m)|≤C_m||w||². Also ∂_m F(p_m)=atanh(m)−g_Nχm and |atanh(m)−m|≥|m|³/3. Therefore

 |∂_m F|+||∇_w F|| ≥ |m|³/3−χΛ N^(-1/2)|m|+(0.08−C_m δ)||w||,

with 0.08−C_m δ>0.0793. The Euclidean tangent gradient norm is bounded below by a fixed q-dependent multiple of this sum, because p=p_m+E y is a fixed invertible orthogonal/soft coordinate chart. This derives (D) **for the displayed rational ε**, without a new unverified third-derivative bound. Outside C_N, the λ term divided by the |m|³ term is at most χΛ/√log N when |m|≥a_N; when |m|<a_N but ||w||≥b_N the fast term controls. The stated N^-3/4(log N)^3/4 lower bound follows, up to fixed constants, for N sufficiently large. Outside C_N the λ term is eventually smaller than a fixed fraction of the RHS's leading term, hence ||∇_T F||≥c' N^(-3/4)(log N)^(3/4). All departure labels have positive occupation, uniformly ≥1/12−ε/12−δ>0, and therefore all q(q−1) elementary directions are available for large N. Select an i of maximal ∂_i F and j of minimal ∂_j F. The discrete step to p+(e_j−e_i)/N reduces F by at least

 c'' N^(-7/4)(log N)^(3/4)−C/N²>0.

By the fixed cylinder/sublevel separation, this step remains in the domain. Repetition must hit C_N, since the count lattice has finitely many points and F strictly decreases. It takes ≤(N+1)^12 steps. Follow the central route to n_*.

Consequently the ENTIRE local type sublevel is connected for all sufficiently large N in 12ℕ, and D_N is precisely the full displayed sublevel. Crucially the descent path has F(noncentral intermediate)≤F(start), while the central route has F≤C log N/N. Since (C) gives F(start)≥−C_Λ/N, every connecting root path x↝n_* obeys the **canonical path energy ceiling**

 max_{z on path(x,n_*)}F(z)≤F(x)+C_*(1+log N)/N. (E)

This stronger-than-connectivity inequality is the input needed by MT-004. The constants are analytic existence constants, not invented numerical spectral gaps or calibrated finite-N thresholds.

## D. Probability coverage and exit layer

For iid initial q12 labels: E m²=1/N, E||w||²≤1/N; Markov's inequality gives

 P_iid(p_0∉C_N)≤N^(-1/2)(log N)^(-1/2)+1/log N=o(1).

Since C_N⊂D_N for N large, μ0_N(D_N^c)=o(1). For any inner boundary departure type n∈D_N with an actual outgoing count edge to D_N^c, neither p nor its neighbour is near the outer cylinder; the neighbour must have F>h. Uniform bounded gradients on C give |F(p')−F(p)|≤C/N. Thus F(p)≥h−C/N≥b=h/2 for all N sufficiently large. This is the paid exit-layer estimate; it does not invoke a global metastable saddle.

**Boundary note:** D_N is a connected sublevel of the exact finite-type F, not a physical wall. Its analytic definition is the same under HB/MET/BAR; only actual jump rates differ.

## Proof obligation and arithmetic scope

All displayed geometric constants ε=1/100, δ=1/10000, h=10^-10 and b=h/2 have explicit rational margins conditional on the inherited *exact* spectral inequalities g_c>5 and γ>0.0186. The Hessian/gradient and connectivity argument is analytic, not a machine-checked global interval proof; N0(Λ) and the canonical-path energy ceiling constant C_* exist but are **not numerically certified**. This suffices for an asymptotic polynomial-mixing existence theorem, not a finite-N confidence claim. If a future review finds a logical error in the discrete gradient/central-corridor construction, MT-004/finite-window must be demoted; MT-002's trapping inequality still follows from any separately admitted positive-boundary domain.
