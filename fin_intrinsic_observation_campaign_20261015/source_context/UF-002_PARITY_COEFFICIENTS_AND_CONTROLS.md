# UF-002 — Microscopic parity characteristics and fast selection (proof)

Let R_{ij}^{(N)}(p,g) be the off-diagonal per-departure-copy rate obtained from h_j=g[A(p-e_i/N)]_j. The **unchanged** rules are HB: softmax_j(h); MET: exp(min(0,h_j-h_i))/12; BAR: (2/12)/(1+exp(-(h_j-h_i))). The actual count rate n_i R_{ij}^{(N)}. In the limit replace p-e_i/N with p. All three limiting rules have strictly positive off-diagonal rates in the interior, are locally Lipschitz in p and g, and the finite-N minus limiting drift/variance is uniformly O(1/N) on compact g (the uniform constant follows from bounded A and bounded derivative/Lipschitz constant of each rate).

For the parity observable m=sigma^T p, each actual jump i->j changes m by (sigma_j-sigma_i)/N (0 or ±2/N). Exactly:

L_N m = sum_{ij} p_i R_{ij}^{(N)}(sigma_j-sigma_i),
N L_N[m,m] = sum_{ij}p_i R_{ij}^{(N)}(sigma_j-sigma_i)^2.

The g=0 per-copy off-diagonal rates are identically 1/12 for all three declared rules, with *no re-tuning*. At p=u, g=g_c, all three have edge intensities c_ij=1/144 for i≠j, parity drift 0 and **N d<m>_t/dt=2+O(1/N)**. Their parity generator jumps, at critical space/time scaling, are ±2 N^(-3/4); total accelerated jump intensity O(N^(3/2)), so absolute third characteristics are O(N^(-3/4)).

## Exact parity line, cubic drift

Since A p_m=chi m sigma, HB's limiting parity drift on p_m is tanh(g chi m)-m. For BAR, a parity flip from source sign s has rate (1-s tanh(g chi m))/2, hence BAR has the exact **same** parity-line limiting drift (not necessarily the same drift off p_m or the same finite-N corrections).

For MET, put t=g chi m. A source with parity s flips at rate (1/2)min(1,exp(-2st)), so

b_m^MET(p_m,g)=(1-m)/2 min(1,e^(2t)) - (1+m)/2 min(1,e^(-2t)).

This is an odd function of m but has a kink in higher derivatives at m=0. At k=g chi close to 1, the m>=0 expansion is (k-1)m+k(1-k)m^2+[(2/3)k^3-k^2]m^3+O(m^4); reflect oddly for m<0. At k=1 the cubic coefficient is -1/3. The O(|k-1| m|m|) contribution is of lower critical order. Do NOT assert C3 regularity of MET rates.

## A uniform selection-rule proof covering the MET kink

For interior p define F_g(p)=sum p_i log(12 p_i) -g p^T A p/2 and c_ij(p)=p_i R_ij(p). Finite q and strictly positive edges give the **exact nonlinear identity**

  log(c_ij/c_ji) = partial_i F_g - partial_j F_g.

For each unordered edge form the positive logarithmic-mean conductance ell(c_ij,c_ji)=(c_ij-c_ji)/(log c_ij-log c_ji), continuously extended to equal values. K(p,g)=sum_{i<j} ell(c_ij,c_ji)(e_j-e_i)(e_j-e_i)^T is positive on the occupation tangent and satisfies **b(p,g)=-K(p,g) grad F_g(p)**. The function K is locally Lipschitz near the uniform interior point even for MET: each edge c is positive locally Lipschitz and log mean is smooth for positive arguments. It is NOT necessary that the MET rates be C3.

Under label shifts by 2, p_m is fixed and both K(p_m,g) and Hess F(p_m,g) commute with the shift representation. On T the only shift-fixed direction is sigma; the ten-dimensional fast representation has no invariant covector. Therefore b(p_m,g) is a multiple of sigma and **sigma^T K(p_m,g) Hess F(p_m,g) w=0** for every w in the fast space.

At p_m, modulo a constant vector,

   grad F_g(p_m) = [atanh(m)-g chi m] sigma,

which is bounded by C(|m|^3+|g-g_c||m|) for |m| small. Expanding grad F to first order in fast w and using local Lipschitz continuity of K (with the above symmetry cancellation) gives the crucial estimate, for all three kinetics,

 |b_m(p_m+w,g)-b_m(p_m,g)| <= C[||w||^2+||w|| (|g-g_c||m|+|m|^3)].

No dangerous O(|m| ||w||) or O(||w||) term remains. This is a local, fixed-model, interior estimate, not a claim of exact finite-N parity Markov closure.

The parity-line drift satisfies uniformly for bounded critical lambda and small |m|:

 b_m(p_m,g)=chi(g-g_c)m-m^3/3+O(|m|^4+|g-g_c|m^2+|g-g_c|^2|m|).

For HB/BAR the O(m^4) bound can be improved, but improvement is unnecessary. With |m|<=R N^-1/4 and g-g_c=lambda N^-1/2, N^3/4 times every displayed remainder tends to zero. Likewise, N^3/4*O(1/N)=O(N^-1/4) pays the finite-N LOO difference without altering kinetics.

## Noise and control scope

For all rules the bracket coefficient B(p,g)=sum p_i R_ij(p)(sigma_j-sigma_i)^2 is continuous near (u,g_c), equals 2 there, and differs from 2 by at most C(|m|+||w||+|g-g_c|). The **microscopic jump bracket**, not an auxiliary Gaussian field, produces sqrt(2) noise. No Metropolis/Barker fold coefficients are imported from the distinct fold g_f. Same leading critical process is an asymptotic local class-level statement for three supplied rules only.
