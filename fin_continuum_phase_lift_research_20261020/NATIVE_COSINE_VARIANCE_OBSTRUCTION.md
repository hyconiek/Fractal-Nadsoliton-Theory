# CPL-005: a positive-mass native moment-cone rejection of the frozen lift

Date: 2026-10-20. This adverse result is retained; no a(h), angle coefficient,
phase threshold or clock is retuned to evade it.

## 1. A stronger necessary native inequality than total trace

Use the ACTUAL integer k3 features C(n),S(n) and whitened coordinate

    eta_C=sqrt(2gamma3) C(n)/sqrt(N).

Every opposite-parity single-copy move changes C by exactly +/-1, whether
upward or downward. Its squared soft jump is 4N^-3/2. Same-parity moves
contribute nonnegative extra cosine variance. Thus under the KR critical
clock, for every state and every nonnegative single-copy generator,

    A_N=Gamma_(sqrt(N)L)(eta_C,eta_C),
    A_N/sqrt(N)>=(gamma3/2)G_N.                       (1)

This is a necessary face of the local moment cone. It requires neither
Gibbs, drift matching, irreducibility nor a numerical approximation.
The actual inherited spectral boxes give gamma3>1/8.

## 2. The frozen CPL coefficients fail (1) on an open set

Keep kappa_N=sqrt(N) in the ALREADY FROZEN generator. In phase A the
normalized cosine variance is

    d_N=2a r_h^2 cos^2(theta)
         +2r^2 sin^2(theta)/(1+r^2)
         +(g_A r_z^2/sqrt(N))cos^2(theta).           (2)

Consider the fixed set U in the existing chart:

    phase A, 19/20<h<1,
    distance(theta,pi Z)<1/20,
    z arbitrary, v arbitrary.

On U, |sin(theta)|<1/20 and 1-h<1/20. The actual radial bound before
discarding its boundary factor is

    2a r_h^2<=[p/(1-p)^2](1-h)^2<=40(1-h)^2/9.

The first two terms of (2) are therefore at most 1/90+1/200=29/1800.
For sqrt(N)>=20 the last term is at most 16/(27*20)=4/135. Consequently

    d_N<=247/5400<1/20,
    (gamma3/2)g_A> (1/8)(1/2)(8/5)=1/10.            (3)

So the CPL cosine coefficient is short by MORE THAN 1/20 on U, despite
its passing total-trace growth. The other eight OU directions cannot
repair this: they are orthogonal to the k3 cosine coordinate.

## 3. Positive stationary mass and an actual L1 obstruction

Under the exact continuum stationary law, U has mass

    mu_A*(1/20)*(1/(10pi))>1/1600,

because mu_A>=1/2, uniform angle gives total interval length 1/5, and pi<4.
The inequalities are strict; its boundary has zero measure.

Let U_N be the same chart-defined set for the actual MT count coordinates.
The admitted static joint limit gives pi_D,N(U_N)->nu(U)>1/1600. This is
a bounded-probability transfer, not an unbounded moment passed through TV.
For ANY candidate count generator, integrate (1) and use (3). With g_P the
same frozen target soft coefficient, its errors must satisfy

    liminf { E_piD |A_N/sqrt(N)-d_N|
                 +(gamma3/2) E_piD |G_N-g_P(Z_N)| }
       >=1/32000>0.                                (4)

The comparison d_N in (4) is its pulled-back field on actual count states;
values off U_N can be defined by the same formula in the two charts.
The zero-mass axis/interface conventions have no effect on the inequality.

This rejects direct L1 local-characteristic transfer of THIS CPL family
to the native single-copy graph. It is stronger than saying such a transfer
has not yet been proved. The result does not reject the original phase
partition or all possible fast dynamics with the desired soft profiles.
It also does not rule out a native SOFT process law without matching these
particular fast characteristics. Full finite-N flow/IBA remains open.

The analysis identifies a concrete new obligation: any future native
construction must satisfy componentwise jump geometry, potentially with
different fast reflection/interface behavior. No replacement is adopted
or automatically commissioned after this adverse result.
