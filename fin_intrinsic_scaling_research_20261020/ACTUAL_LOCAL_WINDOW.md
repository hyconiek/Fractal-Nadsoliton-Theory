# IS-002: actual unaugmented source mixing and local window

Date: 2026-10-20. Conditional on exact model/spectral/MT inputs and the
authorized entropy-source postulate; NOT primary native FIN.

Let P_N be the original count source, B_N the explicit polynomial bound
from IS-001, and P_D its reflection of the ALREADY SELECTED global source.
Every internal edge is retained; no allowed-destination renormalization or
new entropy maximization on pi_D occurs. Detailed balance gives invariant
pi_D=W_N 1_D/W_N(D). Actual MT connectivity applies because all legal source
edges are positive.

## 1. Both ends of the discrete spectrum are paid

For real f, (5) of IS-001 gives

    E_(I-P_D)(f,f)>=B_N^(-1) E_Q_D(f,f)
                    >=a_* B_N^(-1) E_HB_D(f,f),
    1-lambda_2(P_D)>=a_* gamma_*/[B_N(N+1)^K].        (1)

Here a_*=q exp(-2H); gamma_*,K are the actual admitted HB local-gap
constants. The form uses the same nonuniform pi_D throughout. No uniform
measure, native parity drift or constant fast/full gap is assumed.

IS-001 also gives P_D(x,x)>=1/B_N. Since B_N>=d_N>=2, write

    P_D=B_N^(-1)I+(1-B_N^(-1))S_D

with stochastic reversible S_D. The identity
<f,(I+S_D)f>_pi=(1/2)sum_xy pi_x S_D(x,y)(f_x+f_y)^2>=0 pays

    lambda_min(P_D)>=-1+2/B_N.                       (2)

Thus the absolute nonconstant spectral gap is at least

    eta_N=B_N^(-1)min(2,a_* gamma_*/(N+1)^K).        (3)

It is polynomially bounded below. Mere aperiodicity at each finite N, or
the upper side (1) alone, would not be an adequate discrete mixing proof.

## 2. Actual iid density, exit and one-time preparation exception

With nu=mu0(.|D), the exact inherited density calculation gives
||dnu/dpi_D-1||_2<=C(N+1)^6. Therefore

    TV(nu P_D^k,pi_D)<=C(N+1)^6 exp(-k eta_N).        (4)

The actual parent P_N remains unreflected. PSD A7 and eventually g_N>=0
give mu0<=W_N. Detailed balance makes W_N P_N=W_N, hence
mu0 P_N^j<=W_N. Its outward-step probability is <=1, so the actual
departure-event compensator gives

    Pr_iid(tau_D^steps<=k)
      <=delta_N+k (N+1)^q exp(-Nb),
    delta_N=mu0(D^c)->0.                             (5)

This residence identity is reused from SQ; the new conclusion is the
mixing/window gate. Couple reflected nu and actual mu0 paths: on initial
inside states match them until exit; on initial outside states retain the
exception. The coupling failure is at most Pr(tau_D<=k), paying delta once.

Combining (4),(5), for every integer k>=0,

    TV(Law_iid(n_k),pi_D)
      <=delta_N+k(N+1)^q exp(-Nb)
                      +C(N+1)^6 exp(-k eta_N).       (6)

Choose k_loc=ceil((N+1)^(K+4)). Since B_N=O(N^2) and K is a fixed
nonnegative MT exponent, k_loc eta_N grows at least c N^2 eventually.
Thus for every fixed 0<zeta<b,

    sup_(k_loc<=k<=floor(d_N exp(zeta N)))
      TV(Law_iid(n_k),pi_D,N)->0.                    (7)

This is the actual unaugmented hypothesis source on N interacting
coordinates, not a uniform-spectator/Poissonized replacement. No fitted
clock was used. If one writes t=k/d_N, burn-in is polynomial in t as well.
No useful numerical N0/constants or physical duration follows.

## 3. Useful but scoped future-path consequence

At a deterministic late origin k_N in the window, maximally couple its
count law to pi_D, then propagate the same P_N. Whole future-path TV costs
the initial preparation error once, not once per future observation. Over
m further intrinsic steps the reflected stationary comparison adds at most
m pi_D(B_D). For polynomial m this vanishes exponentially under MT.
It compares this source to ITS OWN stationary reflected histories, not to
HB or SQRT histories and not to any yet unpaid critical SDE.

The static conditional Gibbs quartic/Gaussian law remains inherited from MT.
Equation (7) does not transfer unbounded iid moments, posterior/information,
global mixing/escape prefactors, physical energy or a natural source law.
