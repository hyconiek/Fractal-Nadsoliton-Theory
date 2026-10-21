# Oscillatory memory transfer with the causal shape term retained

Date: 2026-10-21. Conditional on the actual NM kernel proof. Same-author
analytical research, not an independently checked full-PDE envelope limit.

## OM-001. Uniform spectral gaps in the frozen family

For nu in the contract's V,

    min(|nu-2|,|nu+2|)>1/5,
    min(|nu-sqrt3|,|nu+sqrt3|)>1/2.                       (1)

The first follows from sqrt5>11/5. For the second, sqrt3<7/4 implies
(sqrt3+1/2)^2<5. Frequencies0 are also separated. These are inherited
source-band facts, not a new mass choice. All constants below are uniform
on V and may depend on the fixed NM kernel and slow horizon S.

## OM-002. Continuum primitive avoids the divergent absolute first moment

Define the ordinary tail

    P_c,nu(t)=integral_t^infinity e^(i nu u)k_c(u)du,
    chi_c(nu)=P_c,nu(0).

NM gives k_c in L1 and its precise t^(-3/2) oscillatory tail. Split the
leading sine into phases nu+2 and nu-2. Integration by parts in each tail,
using (1), bounds their integrals by C(1+t)^(-3/2). The O(u^(-5/2))
remainder has absolute tail O(t^(-3/2)). For t<=1 use ||k_c||_1.
Consequently

    |P_c,nu(t)|<=C(1+t)^(-3/2),
    sup_(nu in V) ||P_c,nu||_1<infinity.                 (2)

This is not the false assertion integral u|k_c(u)|du<infinity.
On a finite interval, P_c'=-e^(i nu t)k_c and integration by parts gives
the EXACT identity

    C_c A(t)=chi_c(nu)A(epsilon^2 t)-P_c,nu(t)A(0)
        -epsilon^2 integral_0^t P_c,nu(u)A'(epsilon^2(t-u))du. (3)

Thus the remainder is <=C epsilon^2||A'||_infinity uniformly even for
t<=S epsilon^-2. Causality at t=0 is exact, not an eternal-response
substitution. Complex-valued envelopes, including slow phase changes,
are permitted; no derivative bound for an actual physical envelope is
inferred from this test-function statement.

## OM-003. The shape tail is bounded but not integrable

Put lambda=3, eta=sqrt(lambda), D=lambda-nu^2 and chi_s=d^2/D. Define
bounded algebraic/Abel primitives, NOT ordinary integrals to infinity:

    P_s,nu(t)=chi_s e^(i nu t)[cos(eta t)-i nu sin(eta t)/eta],

    Q_s,nu(t)=-d^2 e^(i nu t)/D^2
                 [2i nu cos(eta t)+(lambda+nu^2)sin(eta t)/eta]. (4)

Direct differentiation proves P_s'=-e^(i nu t)k_s and Q_s'=-P_s.
Also P_s(0)=chi_s. The gap(1) bounds P_s and Q_s uniformly on V.
Neither is assumed decaying or L1.

The same finite-interval identity as(3) holds. Integrate its last term
once more, using Q_s, to get

    C_s A=chi_s A(epsilon^2 t)-P_s(t)A(0)
          -epsilon^2 Q_s(0)A'(epsilon^2 t)
          +epsilon^2 Q_s(t)A'(0)
          +epsilon^4 integral_0^t Q_s(u)A''(epsilon^2(t-u))du. (5)

On t<=S epsilon^-2 this remainder is bounded by
C epsilon^2[||A'||+S||A''||]. The nondecaying initial term P_s(t)A(0)
is retained. For a constant envelope (5) is exact with zero error; at
nu=0 it is d^2 A(0)(1-cos(sqrt3 t))/3. Discarding that oscillation would
give an order-one error, not a small transient.

## OM-004. Combined response, derivative and scaled-drive bound

Let chi=chi_c+chi_s and P=P_c+P_s. Then

    C_(nu,epsilon)A(t)=chi(nu)A(epsilon^2 t)-P_nu(t)A(0)+R(t),

    sup_(0<=t<=S epsilon^-2)(|R(t)|+|R_t(t)|)
      <=C epsilon^2[||A'||+(1+S)||A''||].                (6)

The value bound is(3),(5). For the derivative, use the exact remainder
R=-epsilon^2 integral_0^t P(u)A'(epsilon^2(t-u))du:

    R_t=-epsilon^2 P(t)A'(0)
         -epsilon^4 integral_0^t P(u)A''(epsilon^2(t-u))du.

P_c is bounded/L1 and P_s bounded, so its last absolute integral is
<=C[1+t]||A''||. This proves the original-t derivative bound in(6).
No third derivative, finite absolute first moment or small-t cutoff is
hidden in it.

For a drive r(t)=epsilon^2 e^(-i nu t)A(epsilon^2 t), multiply(6) by
epsilon^2 e^(-i nu t). Since |nu|<=sqrt7, the projected response and its
original-time derivative have O(epsilon^4) error after the displayed
susceptibility AND initial slip are retained. The result concerns
<F,eta>, not a bounded global L2 outgoing-resolvent field.

The constants use the old spectral response: chi(0)=1/70 and, for
nu=Omega>0, Im chi(Omega)=Gamma_old/9. Those are reused susceptibility
identities, not a newly proved averaged mode-decay coefficient or noise.

## OM-C01. Actual continuum-threshold counterexample

Let C0=1/(16sqrt(pi)), the actual NM tail coefficient. At nu=2 the
product e^(i2t)k_c(t) has a nonoscillatory leading term
(C0/2)e^(-i pi/4)t^(-3/2). Therefore

    P_c,2(t)=C0 e^(-i pi/4)t^(-1/2)+O(t^(-3/2)).         (7)

It is not L1. Take A(s)=s. Formula(3) remains an exact finite-interval
identity, but at t=S epsilon^-2 its remainder is

    -2C0 e^(-i pi/4)sqrt(S)epsilon+O(epsilon^2),          (8)

not O(epsilon^2). The bounded nonresonant shape error there is only
O(epsilon^2) and cannot cancel (8). Thus removing the threshold gap from
(6) is FALSE for the actual kernel, not only for an artificial tail model.
Frequency2 lies outside the frozen source band; no a0 is changed to make
this into an adverse claim about the admitted family.

At the shape resonance nu=eta, a constant envelope instead gives
C_s(t)=i d^2 t/(2eta)+a bounded oscillation. There is no finite chi_s or
bounded primitive(4). That resonance is likewise outside the frozen band.

## OM-C02. Uniform C1 envelope control alone is insufficient

Even with base carrier nu=0 IN V, take

    A_epsilon(s)=epsilon^2 exp(-i eta s/epsilon^2).

The envelope and its first slow derivative are uniformly bounded, but
||A_epsilon''||=3/epsilon^2. Its physical-time drive excites the actual
shape resonance. Direct integration gives its shape convolution

    epsilon^2 i d^2/(2eta)
                  [t exp(-i eta t)-sin(eta t)/eta].

At t=S epsilon^-2 the first term has modulus d^2 S/(2eta). The proposed
susceptibility/initial-slip terms are only O(epsilon^2), and the continuum
response to this small input is O(epsilon^2) by ||k_c||_1. The remainder
therefore stays of order one, not O(epsilon^2). This is not a new source
packet: it is a test-family counterexample to dropping the uniform second-
derivative norm in the operator statement. It does not assert that this
family is the actual RP envelope.
The real nonnegative family1+epsilon^2 cos(eta s/epsilon^2) gives the
analogous failure of the uniform sup bound at its resonance peaks; its
constant part is handled exactly by the retained initial slip.

## Exact scope of the advance

This pays the projected memory-operator component of NM's long-time gate.
It does NOT prove that the full nonlinear solution has a slow envelope
with uniform C2 norms, extend its O(epsilon^5(1+t)^4) prediction bound to
epsilon^-2, remove the coherent mode by a nonlinear normal form or prove
an averaged amplitude equation. Full PDE remainders/harmonics/Sobolev
control remain open. No new microscopic law, noise, initial packet,
physical time, F01/observer/QM/GR/ToE or selector-role claim is supplied.
