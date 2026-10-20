# CPL-002: existence, stationary Gibbs law and actual nonaveraging

Date: 2026-10-20. Conditional analytical construction; not native counts.

## 1. The soft process exists for all comparison times

The inherited profiles are smooth, even and positive. Their tail representation
and monotonic p_A,p_B give

    8/5<=g_A(z)<2<g_B(z)<=8/3 for every finite z,
    J_A+J_B=rho, p g_A+(1-p)g_B=2.                    (1)

For z>=0 use J_P=integral_z^infinity (u^3/3)q_P(u)du and
integral_z^infinity (u^3/3)rho(u)du=rho(z). The strict inequalities follow
from strict monotonicity of p on every positive-length tail; evenness covers
negative z. These are reused/verified FP profiles, not a new cap theorem.

The Z equation has locally Lipschitz coefficients and bounded nonzero
diffusivity. Its generator on z^2 is -2z^4/3+g_P(z)<=8/3. Stopped Ito
estimates therefore exclude finite explosion. It has a unique global
strong solution in each chart.

## 2. No hidden endpoint loss in H

The H coefficients are smooth inside (0,1). Its scale derivative is 1/a;
one primitive is

    s(h)=log h-log(1-h)+1/(1-h), s'(h)=1/[h(1-h)^2]. (2)

It diverges at both endpoints, and partial_h(a s')=0. Stopping s(H) on
compact subintervals yields the ordinary scale hitting-probability formula.
Letting the lower endpoint tend to 0, or the upper one to 1, makes the
probability of hitting that endpoint before a fixed interior return zero.
A countable rational-time argument covers a putative finite endpoint hit:
each continuous hitting path has a rational time after its last return to
a fixed interior level. Thus neither endpoint is reached from the interior.
Local strong existence/uniqueness extends for all finite times. This proof
does not merely assert that a formally invariant density prevents explosion.

On a finite interval H has a positive distance from both endpoints almost
surely. Z has a finite maximum; hence r_P(Z,H) stays finite and positive,
the angle equation is well-defined, and v is the usual nonexplosive OU.
The pulled-back state paths are continuous on every finite horizon.

## 3. Pay invariance and reversibility, not just a formal PDE

Since (g_P q_P/2)'=J_P'=b q_P, integration by parts gives the symmetric
form in chart P with measure (2):

    E_P(F,F)=integral [(g_P/2)F_z^2
          +kappa(a F_h^2+c_P F_theta^2+||grad_v F||^2)] dnu_P. (3)

For Z and H, one direct conservativity/domain payment is one-dimensional
compact exhaustion. Reflect Z on [-R,R] and H on [delta,1-delta]. Their
smooth divergence-form Neumann operators are conservative and self-adjoint
in the truncated densities q_P and uniform dh. The true solutions couple
until the first artificial boundary hit. Nonexplosion and endpoint
inaccessibility make that probability tend to zero on each fixed horizon;
truncated initial laws converge in TV. Thus stationarity and reversibility
pass to Z,H. The independent OU coordinates are reversible by their centred
Gaussian covariance exp(-kappa|s-t|).

To include the VARIABLE angle clock without an unsupported product claim,
condition on the entire reversible base path (Z,H,v). Angle transitions
over [s,t] are circular heat kernels with time parameter
kappa integral_s^t c_P(Z_r,H_r)dr. Each kernel is symmetric and preserves
uniform angle. Under reversal of the base path these same clock intervals
reverse order. The conditional finite-cylinder density, uniform initial
angle times the product of the symmetric heat kernels, is unchanged by
reversing all angles. Hence the full process is stationary and reversible.
This also supplies the conservative symmetric realization of (3).

The angle coefficient depends on Z,H but not theta, so it creates no extra
stationary flux or derivative in those other coordinates. Mixture of the
two closed reversible charts and their bijective pullback proves exact
reversibility with the FULL rho times fast-Gaussian law. No hidden extra
degree of freedom, physical wall or natural mediator is inferred.

## 4. Exact failure of the continuum averaging implication

The generator on an observable f(z) is b f'+(g_P/2)f''. Thus its actual
soft quadratic-variation density is G(z,Y)=g_(P(z,Y))(z), whereas

    E[G|Z=z]=2 EXACTLY.                              (4)

Let M=E_rho[Z^2 exp(-Z^2)]>0. Quartic-tail integration of
(z exp(-z^2)rho)' gives
E[(Z^4/3-1)exp(-Z^2)]=-2M. Integration by parts of J_P then gives

    integral q_A(g_A-2)=-M/2,
    integral q_B(g_B-2)=+M/2.                         (5)

Every path stays in its intrinsically defined phase. By (1), G-2 has one
strict sign along the WHOLE path. Conditional stationarity and boundedness
therefore give the exact identity, for every S>0,

    E sup_(t<=S)|integral_0^t(G(Z_s,Y_s)-2)ds|=S M>0. (6)

This is not just a positive sufficient Poisson cost. It realizes the old
IF-feasible obstruction as an ACTUAL CONTINUUM process, not as a finite-N
family. Its law of Z and (6) are unchanged for every kappa>0, even if a
sequence of comparison kappa values grows.

For a coarse fully analytic bound, the numerator defining M, restricted to
1/2<=|z|<=1, exceeds (1/4)exp(-13/12)>1/16. Here e<3 and
exp(1/12)<12/11 give exp(13/12)<4. The quartic normalizer is less than
6+2integral_3^infinity exp(-z)dz<8. Hence M>1/128. This is a comparison-
model lower bound, NOT a finite-N IBA refutation or a measured quantity.
