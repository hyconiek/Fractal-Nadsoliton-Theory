# VD: exact polynomial virial coercivity with one threshold defect

Date: 2026-10-21. Conditional analytical proof, same-author audit.
A continuum block and actual critical-time estimate, not a full modal FGR law.

## 1. Fixed weight, source-compatible phase space

Use the unchanged RP real/parity preparation and parameter band.
QP phase positions/momenta are (a_i,p_i); their position derivatives
include QP's forcing and are not silently identified with p_i.
Remove the actual bound profiles with Pc. Use the normalized Darboux map

    T=partial_x^2-3K partial_x+3K^2-1,
    P=(-partial_x^2+1)(-partial_x^2+4),
    Wcal=T P^-1/2, Wcal*=P^-1/2(partial_x^2+3K partial_x+2),
    w_i=Wcal* Pc a_i, pi_i=Wcal* Pc p_i.

The paid RP resolution gives Wcal* Wcal=1, Wcal Wcal*=Pc, and intertwines
H/L with the FREE operators -partial_x^2+m_i^2,
m_phi=2, m_psi=sqrt(4+w). Both maps preserve realness and parity.

Fix ONCE

    varphi(x)=arctan x, zeta=(1+x^2)^(-1/2), W=zeta^3=<x>^-3,
    q(w)=integral_(-1)^1 zeta(x)w(x)dx.                  (1)

This polynomial choice matches QP's proved polynomial moments. No weight
scan, mass fit or physical premise is introduced. q is a bounded compact
average, not a pointwise momentum trace, an added state or a new law.

The Darboux maps preserve the weighted phase norms used below.
Here is the domain proof. Let R_j=(-partial_x^2+j^2)^-1/2.
Each multiplier R_j and partial_x R_j has bounded derivatives of all
orders. Commuting integer x moments in Fourier proves boundedness on
positive-polynomial weighted L2/H1. Duality gives negative weighted L2;
commutation with partial_x gives negative weighted H1. Multiplication by
K and its bounded derivatives preserves these norms. Explicitly,

    Wcal=(partial_x R1)(partial_x R2)
          -3K(partial_x R1)R2+(3K^2-1)R1R2,
    Wcal*=(partial_x R1)(partial_x R2)
          +3(partial_x R1)R2 K-3R1R2 K'+2R1R2.

Thus both are bounded for weights <x>^+/-N at H1/L2 regularity.
On Pc, inverse identities give norm equivalence. No T* derivative loss
or unsupported weighted singular-integral theorem is invoked.

## 2. Exact coercivity of the free virial form

For real w in H1 define

    B(w)=integral varphi' w_x^2-(1/4)varphi''' w^2.

Set u=zeta w. Integration by parts gives exactly

    B(w)=||u_x||2^2+(1/2)integral [(x^2-1)/(1+x^2)^2]u^2. (2)

Only |x|<1 contributes a negative potential, whose magnitude is <=1/2.
Write q=integral_(-1)^1 u and ubar=q/2. An elementary mean inequality is

    integral_(-1)^1 |u-ubar|^2
         <=(1/2)integral_(-1)^1(1-s^2)|u_x(s)|^2 ds
         <=(1/2)||u_x||2^2.                             (3)

For completeness, use
integral |u-ubar|^2=(1/2)integral_(x<y) |u(y)-u(x)|^2 dxdy,
then Cauchy-Schwarz on [x,y] and Fubini. The coefficient of |u_x(s)|^2 is
(1/2)integral_(-1)^s integral_s^1(y-x)dydx=(1-s^2)/2.
No numerical eigenvalue or scalar-paper coercivity premise is imported.

Equations(2)--(3) imply

    B(w)>=(3/4)||u_x||2^2-q^2/4.                        (4)

To control local position energy, average u(x)-u(y) over -1<=y<=1:

    |u(x)|^2<=q^2/2+2(|x|+1)||u_x||2^2.

Using integral(1+x^2)^-2=pi/2 and integral|x|(1+x^2)^-2=1 yields
||zeta^2 u||2^2<=(2+pi)||u_x||2^2+pi q^2/4.
Since |(zeta^2)'/zeta^2|<=1,

    ||Ww||H1^2<= (8+3pi)||u_x||2^2+(3pi/4)q^2
                 <=(122/7)||u_x||2^2+(33/14)q^2.        (5)

The last rational bound uses the already paid elementary pi<22/7.
Combining(4)--(5) gives the explicit continuum coercivity

    B(w)+(1/2)q(w)^2 >=(1/24)||Ww||H1^2.                (6)

Indeed the two residual coefficients are1/42 and17/112, both positive.
For odd w, q=0. For even w it is a REAL defect, not a parameter to set0.
This is an instantaneous inequality, not yet integrated control of q.

## 3. Forced phase virial including momentum

Let a free phase solve

    w_t=pi+f1, pi_t=w_xx-m^2 w+f2, 2<=m<=sqrt(23/4).

For smooth solutions define

    I=integral (varphi w_x+(1/2)varphi' w) pi,
    J=integral W^2 w pi, delta=1/816, V=I-delta J.

Integration by parts gives EXACTLY

    I'=-B(w)+E_I,
    E_I=integral (varphi f1_x+(1/2)varphi' f1)pi
                       +(varphi w_x+(1/2)varphi'w) f2,
    J'=||Wpi||2^2-||Ww_x||2^2
          +(1/2)integral(W^2)''w^2-m^2||Ww||2^2+E_J,
    E_J=integral W^2(f1 pi+w f2).                       (7)

No derivative of f2 or pointwise momentum trace is used.
For a=W^2=(1+x^2)^-3,

    a''/a=(-6+42x^2)/(1+x^2)^2, |a''|<=12a,
    |W'/W|<=3/2.

Consequently J'>=||Wpi||2^2-17||Ww||H1^2+E_J:
the largest L2 coefficient is23/4+6+9/2=65/4<17.
Since1/24-17delta=1/48, (6)--(7) yield

    -V'+(1/2)q^2
      >=(1/816)[||Ww||H1^2+||Wpi||2^2]
                             -|E_I|-delta|E_J|.        (8)

I,J are finite and bounded by C(||w||H1+||pi||2)^2.
The identities first hold on a smooth core and extend by energy
approximation for f1 in H1, f2 in L2. q is not differentiated.

## 4. Apply it to the ACTUAL coupled source

QP's exact phase system, after Pc and Wcal*, has precisely the form
in section3, for i=phi,psi, with actual forcing split into loc and ge3.
Its source/coordinate bounds give uniformly in every finite model time

    ||(w_i,pi_i)||_(H1 x L2)<=C epsilon,
    ||(<x>^3 f1_loc,<x>^3 f2_loc)||_(H1 x L2)<=C epsilon^2,
    ||(f1_ge3,f2_ge3)||_(H1 x L2)<=C epsilon^3.           (9)

The positive moments follow from QP and the weighted Darboux bounds.
The phase conversion uses f1=(R^++R^-)/2 and
f2=Omega(R^+-R^-)/(2i), then Wcal* Pc. Omega maps positive weighted H1
to weighted L2 by its bounded symbol-derivative commutators.
The original noncanonical position-velocity distinction is retained.

For the local forcing,
|E_I,loc|+delta|E_J,loc|<=C epsilon^2 E_i(t), where
E_i^2=||Ww_i||H1^2+||Wpi_i||2^2. Pair a negative weight on the field
with the derived positive weight on the source. All varphi coefficients
are bounded and |W'/W|<=3/2. Young's inequality absorbs half of(8)'s
energy term and leaves C epsilon^4. For cubic forcing, use GLOBAL
C epsilon energy times C epsilon^3 source; its contribution is also
bounded by C epsilon^4. No positive moment of cubic radiation is assumed.

Oddness gives q_phi=0. Summing (8), and writing
q_psi(t)=q(w_psi(t)), integration proves the ACTUAL bound

    integral_0^T [E_phi(t)^2+E_psi(t)^2]dt
      <=C epsilon^2+C epsilon^4 T+816 integral_0^T q_psi(t)^2dt. (10)

The boundary functional is C epsilon^2 by actual global energy.
This estimate is valid for every finite T. In particular, for each fixed S,

    integral_0^(S epsilon^-2) E_rad(t)^2dt
       <=C_S epsilon^2+816 integral_0^(S epsilon^-2)q_psi(t)^2dt. (11)

It is a CRITICAL-TIME estimate with one explicit scalar trace defect,
not a proved small trace or a completed nonlinear damping theorem.

Darboux weighted equivalence and the O(epsilon^2) inverse QP correction
transfer the same form of estimate to the OLD PHYSICAL continuous fields
(Pc eta,Pc eta_t,Pc psi,Pc psi_t), with W=<x>^-3. Its constants, including
the trace coefficient, become source-dependent C after norm comparison.
The coordinate correction contributes at most C epsilon^4 T.
Bound modes themselves are not included in this continuum estimate.

## 5. The trace is the already identified threshold, not a new degree of freedom

Let r0=(3K^2-1)/2. Direct computation gives T* r0=2 and hence

    Wcal* r0=1

as a generalized function. The massive inverse multipliers act on
constants by1/j; the same result follows by weighted approximation.
Thus CD's actual local threshold history Theta=r0(theta1,theta2)
becomes a constant phase vector in the free Darboux frame.

The previously paid deflated actual remainder, together with negative-
polynomial weighted Darboux boundedness, gives

    q_psi(t)=2arsinh(1) theta1(t)+r_q(t),
    |r_q(t)|<=C(epsilon^2+epsilon^3 t).                 (12)

This is a relation to the ACTUAL causal history, not an autonomous ODE
for q or a primitive physical memory. It uses only a compact test function
after the weighted operator; no unjustified global subtraction of a
non-L2 threshold state is made.

## 6. Exact remaining joint problem

All continuum coercivity and forcing terms in (11) are now paid, with
the displayed trace left visible. The bound q<=C epsilon alone permits
an order-one integral on the critical interval and does not close (11)
at order epsilon^2. A genuinely sharper trace/history budget is needed.

Furthermore the actual corrected-mode action increments still need their
own error identity before the positive frozen three-channel work form
can become a damping law. This proof does not assert that positivity of
those coefficients controls every mixed modal term. No full joint virial,
o(epsilon^2) action-law error, natural action/source, clock or ToE closure.

