# VA: virial, source budget and weighted absorption audit

Date: 2026-11-07. Parent: b5e77b1dcbd99bcc0827a759c234718f164f93c1.
Same-author analytical audit; independent review remains outstanding.
Disposition: ACCEPT_AT_STATED_DEPENDENCY_BOUNDARY.

## 1. Scope and inherited inputs

We rederive VD sections 2--4 and LF's cumulative-to-weighted argument,
then check the ML/CS substitution. The supplied model, preparation and
w in [5/4,7/4] are unchanged. E is the physical continuous H1-position /
L2-velocity norm with weight W=<x>^-3, excluding bound modes.

Inherited inputs are SP global smallness, the Darboux resolution and its
weighted equivalence, QP's exact coordinate/source bounds, the localized
response estimates used in LF/ML, and the subsequent TP2 cubic trace
estimate. This audit checks how those bounds enter the virial argument;
it does not independently reprove all of those input theorems. No action
law or persistent O(epsilon^2) shape-mode hierarchy is assumed.

## 2. Coercivity and the real trace defect

Set phi=arctan x, zeta=<x>^-1, u=zeta w and
q=integral_(-1)^1 u. For the scalar free virial form,

    B(w)=integral phi' w_x^2-(phi'''/4)w^2
        =integral u_x^2+[(x^2-1)/(2(1+x^2)^2)]u^2.     (VA1)

The equality is integration by parts, not pointwise equality of densities.
The discarded cross term is an explicit total derivative and is checked
symbolically in the companion script.

On [-1,1] the mean is q/2. For an interval of length two,

    integral |u-q/2|^2=(1/2)integral_(x<y)|u(y)-u(x)|^2 dxdy.

Cauchy-Schwarz on [x,y] and Fubini produce the coefficient
(1/2) integral_(-1)^s integral_s^1(y-x)dydx=(1-s^2)/2.
Consequently the local variance is at most ||u_x||2^2/2, and

    B(w)>=3||u_x||2^2/4-q^2/4.                       (VA2)

Averaging u(x)-u(y) over [-1,1] gives
|u(x)|^2<=q^2/2+2(|x|+1)||u_x||2^2. Integrating with <x>^-4
and using |(zeta^2)'/zeta^2|<=1 yields

    ||Ww||H1^2<=(8+3pi)||u_x||2^2+(3pi/4)q^2
               <=(122/7)||u_x||2^2+(33/14)q^2.

Thus B(w)+q^2/2 >= ||Ww||H1^2/24. The two remaining positive
coefficients are 1/42 and 17/112. Odd parity sets q_phi=0; even parity
does not set q_psi=0. No coercivity without that defect is exported.

## 3. Forced phase identities and spatial boundary terms

Let w_t=p+f1 and p_t=w_xx-m^2w+f2. Here p is a phase momentum;
it need not equal w_t. Define

    I=integral (phi w_x+phi' w/2)p,
    J=integral a w p, a=W^2, V=I-J/816.

The homogeneous identities are I'=-B and
J'=||Wp||2^2-||Ww_x||2^2+integral a''w^2/2-m^2||Ww||2^2.
For compact smooth fields the exact spatial fluxes are

    flux_I=phi p^2/2+phi w_x^2/2+phi' w w_x/2
                         -phi'' w^2/4-m^2 phi w^2/2,
    flux_J=a w w_x-a' w^2/2.

These identities are checked as differential polynomial identities.
For energy data, cut off and approximate in H1 x L2. Bounded phi and
its derivatives, and bounded a derivatives, make the boundary errors
vanish in the integrated identity. No pointwise spatial decay assertion
or deletion of an arbitrary finite-box flux is needed.

The exact forcing terms are

    E_I=integral (phi f1_x+phi' f1/2)p+(phi w_x+phi' w/2)f2,
    E_J=integral a(f1 p+w f2).

Both f1 and f2 are retained. No derivative of f2 or differentiation of q
is required. The identities extend on finite time intervals for the
recorded H1 x L2 source class by energy approximation.

For a=<x>^-6, a''/a=(-6+42x^2)/(1+x^2)^2. Its range is
[-6,147/16], so the historical bound |a''|<=12a is valid. Also
|W'/W|<=3/2. Writing W w_x=(Ww)'-(W'/W)Ww gives

    J'>=||Wp||2^2-17||Ww||H1^2+E_J,

because m^2<=23/4 and 23/4+6+9/2=65/4<17. Therefore

    -V'+q^2/2 >= (||Ww||H1^2+||Wp||2^2)/816
                                      -|E_I|-|E_J|/816. (VA3)

The position coefficient before weakening is 1/24-17/816=1/48.
The negative mass term in J' is explicitly included. The mass terms in
I' cancel as a total derivative; they are not omitted from the PDE.

## 4. Complete source and temporal endpoint budget

In the QP/Darboux continuous phase, the inherited bounds are global
phase O(epsilon), positively weighted localized forcing O(epsilon^2),
and unweighted higher forcing O(epsilon^3).

| Contribution | Bound used in VA3 | After time integration |
| --- | --- | --- |
| Initial and terminal V | |V(0)|+|V(T)| <= C epsilon^2 | C epsilon^2 |
| Localized forcing | C epsilon^2 E <= E^2/1632+C epsilon^4 | Half the coercive term plus C epsilon^4 T |
| Higher forcing | Global phase times global source <= C epsilon^4 | C epsilon^4 T |
| Even trace | q_psi^2/2 | Retained explicitly |
| QP-to-physical correction | Global phase correction O(epsilon^2) | Squared local correction C epsilon^4 T |

The last row is a change of observed coordinates, not an additional
forcing term. Weighted Darboux equivalence and a square-of-sum inequality
change constants but no epsilon exponents. The actual QP initial shift
is covered by the endpoint estimate and is not reset to zero.

After absorbing the localized forcing and summing both fields, the
QP-frame trace coefficient is 816. Norm comparison and correction may
increase it to a source-dependent D in physical coordinates. Hence

    F(t):=integral_0^t E(s)^2ds
       <= A epsilon^2+B epsilon^4 t+D Q(t),
    Q(t)=integral_0^t |q_psi(s)|^2ds.                 (VA4)

The defect q is the actual QP compact observation throughout this step.
Its later physical-response decomposition includes the chart correction.
On t=S epsilon^-2, epsilon^4 t=S epsilon^2, not O(1). No unweighted
cubic source is assigned an unproved spatial moment. The estimates do
not require the individual source-response pieces to have global
small energy; that bound is used only for the whole actual phase.

## 5. Weighted cumulative inequality without differentiating VA4

Set beta=2mu epsilon^2, mu>=1. Since F(0)=Q(0)=0,

    integral_0^T exp(-beta t)dF(t)
      =exp(-beta T)F(T)+beta integral_0^T exp(-beta t)F(t)dt.

Substitution of VA4 is legitimate because the coefficients on its right
are nonnegative. The constant part becomes A epsilon^2 exactly. The
linear part becomes B epsilon^4(1-exp(-beta T))/beta. The Q part becomes
D integral_0^T exp(-beta t)|q_psi(t)|^2dt. Thus

    ||E||_alpha,T^2 <= [A+B/(2mu)]epsilon^2
                                      +D||q_psi||_alpha,T^2,
    alpha=mu epsilon^2.                              (VA5)

A cumulative inequality cannot generally be differentiated: for example
F(t)=t+sin(t)/2 is nondecreasing and F(t)<=G(t)=t+1/2 for t>=0,
but F'(0)=3/2>G'(0)=1. This exact negative control is in the script.
Keeping the temporal endpoint in the Stieltjes identity avoids that error.

## 6. Absorption, fixed-S constants and the CS substitution

ML's exact decomposition is q_psi=p_new+q_fb+q_ge3, where
q_fb=q_L+q_2loc and

    |q_fb(t)|<=C epsilon integral_0^t (1+t-s)^-1/2 E(s)ds,
    |p_new(t)|<=C[epsilon^2+epsilon^3 sqrt(1+t)].

The sum defining q_fb is estimated before squaring. Exponential weighting
and Young give ||q_fb||_alpha,T^2<=C_fb||E||_alpha,T^2/mu,
since integral_0^infinity exp(-alpha r)r^-1/2 dr=sqrt(pi/alpha).
The source p_new has squared weighted integral O(epsilon^2): its terms
scale as epsilon^4/beta, epsilon^6/beta and epsilon^6/beta^2.

Insert the square-of-three-summands bound in VA5 and choose once
mu>=max(1,6D C_fb). The feedback coefficient is then at most 1/2.
Absorption gives

    ||E||_alpha,T^2 <= C epsilon^2+6D||q_ge3||_alpha,T^2.

Unweighting at T=S epsilon^-2 costs exp(2mu S). This constant is explicit
in its dependence on S, although D and C_fb are inherited source
constants, not evaluated numerical values. No uniform all-time bound or
physical exponential damping is inferred from the analytical weight.

CS writes q_ge3=q_3+q_hi. Its global quartic forcing estimate gives
integral |q_hi|^2 <= C_S epsilon^3, which is absorbed in C_S epsilon^2
for epsilon<=1. This leaves exactly CS4. TP2's separate structured cubic
trace bound then yields integral E^2 <= C_S epsilon^2. Neither step uses
AW's resulting modal law, so this route is not circular.

## 7. Verdict and next obligation

No blocking defect was found in the audited virial and absorption steps.
The exact constants, both phase sources, temporal endpoints, chart shift
and critical-time powers are accounted for. This pays the upstream audit
requested by AA at the stated inherited-input boundary.

Independent proof review remains open. The next focused analytical audit
should examine LE/QH's even-channel propagator tails and modulated causal
susceptibility, including the threshold term, source regularity and
switch-on transient used in RW4--RW6. Those response estimates were inputs
here, not consequences of VA. SP existence, CM chart construction and CD
coefficient calculations also remain explicitly inherited. This result
adds no natural source, physical units or empirical FIN validation.
