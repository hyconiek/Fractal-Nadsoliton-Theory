# EF: the cubic even-response budget for the finite-packet comparison

Date: 2026-11-05. Base: a12ea9ab9d32dc21fc3049f9a6440663cccb0b37.
Conditional packet-comparison result. Not an approximation theorem for the
original nonlinear RP preparation; not formal verification or peer review.

## 1. Keep the full cubic comparison hierarchy

Use TP's real even packet P_R and PC's odd continuous response Y_R on
0<=t<=T=S R^2, R>=1. In addition retain the generated shape coefficient:

    c_R''+3c_R=-6<g,K P_R^2>,  (c_R,c_R')(0)=0.

The full odd second-order response to this comparison source is
Y_R+c_R g. Define three zero-data EVEN CONTINUOUS responses under L:

    Z_C,tt+L Z_C=-12 P_e(K P_R Y_R),
    Z_S,tt+L Z_S=-12 P_e(K P_R c_R g),
    Z_D,tt+L Z_D=-P_e(P_R^3).                         (EF1)

The cubic comparison response for a fixed signed coefficient g4 is
Z_3=Z_C+Z_S+g4 Z_D. No coefficient is fitted. This is the third-order
source hierarchy for packet data: eta begins at order two, so eta^2 psi
starts at order five; g6 psi^5 also starts at order five. The actual RP
initial data have no such packet, and no long-time nonlinear remainder
for packet data is asserted here.

For q_j=q(Wcal*Z_j), we prove

    integral_0^(S R^2)|q_C|^2dt<=C_S R^4,
    integral_0^(S R^2)|q_S|^2dt<=C_S R,
    integral_0^(S R^2)|q_D|^2dt<=C_S R^3.             (EF2)

At packet amplitude epsilon, all three cubic responses scale as epsilon^3.
For R=epsilon^-1, therefore,

    integral_0^(S epsilon^-2)|q(epsilon^3 Z_3)|^2dt
                                                   <=C_(S,g4) epsilon^2. (EF3)

This meets the desired order for THIS COMPARISON FAMILY. It does not close
the actual nonlinear continuum gate or any modal damping law.

## 2. Dual odd local energy gives a global response bound

PC paid the odd continuous homogeneous estimate

    ||U_H,c(.) X||L2(R;X_-6)<=C||X||X_0.

On the odd continuous sector use the conserved positive energy inner
product, equivalent to X_0. For any continuous phase forcing F=(f1,f2)
and any test phase X, integration by parts in space gives

    |<F,U_H,c(s)X>_energy|
                    <=C||F||X_6 ||U_H,c(s)X||X_-6.

The first-coordinate pairing is <H f1,a>, defined by its H1 quadratic
form; it uses only first spatial derivatives and the bounded potential.
Derivatives of the polynomial weights are bounded relative to those
weights. Thus no H2 source assumption or pointwise velocity trace enters.
Cauchy-Schwarz in time and PC's estimate imply by duality

    ||integral_0^t U_H,c(t-s)F(s)ds||X_0
                                  <=C||F||L2(0,t;X_6). (EF4)

The energy group is unitary in its energy inner product. The inequality
is for the full retarded integral at a fixed t, obtained by testing its
Duhamel integral; no unsupported retarded spacetime theorem is imported.

Insert PC's exact decomposition with Q_c and its initial transient.
There ||Q_c(t)||X_0<=C_S R^-1/2 and ||F_c(t)||X_6<=C_S R^-1.
Consequently

    ||(Y_R,Y_R,t)(t)||X_0
                       <=C_S[R^-1/2+R^-1 sqrt(t)]<=C_S,
                       0<=t<=S R^2.                 (EF5)

Both Q_c(0) and the causal history are retained. The upper bound allows
global radiated energy to be order one for the unit packet comparison;
it does not claim spatial localization of the outgoing response.

## 3. The broad continuum feedback needs no positive spatial moment

PC/TP give ||P_R||infinity<=C_S R^-1/2. With EF5,

    ||-12P_e(K P_R Y_R)||2<=C_S R^-1/2.               (EF6)

The product is allowed to be broad. P_e is bounded on L2 and Wcal* is
unitary on the continuous subspace. The forcing thus enters the free
Darboux phase as a global L2 momentum source of this size. CT's actual
finite-time observation estimate gives

    ||q_C||L2(0,T)
       <=C(1+T)^1/4 integral_0^T ||source_C(s)||2 ds
       <=C_S(1+T)^1/4 T R^-1/2.

Squaring and using T=S R^2 proves the R^4 bound in EF2. There is no
assumption of positive moments or removal of an outgoing wave. The
initial normal-form transient is already included through EF5.

## 4. The shape response is retained and smaller in this budget

TP's exact-packet source expansion and DT's d_S=<g,K r0^2> yield

    -6<g,K P_R^2>
       =-(3/(2R))d_S[A(t/R^2)^2 exp(-2imt)
                         +2|A(t/R^2)|^2
                         +conj(A(t/R^2))^2 exp(2imt)]
                                      +e_c(t),
    |e_c(t)|<=C_S R^-3.                               (EF7)

This scalar pairing is legitimate with TP's negative-weight error because
g has every positive spatial moment. It does NOT promote that local error
to a positive-weight forcing bound for an extended PDE source.

For the shape oscillator lambda=sqrt(3), the carriers 0,+/-2m are all
uniformly separated from +/-lambda on the mass band. The constant-
amplitude modulated sine/cosine primitives are therefore uniformly
bounded. Integrating against the C1 envelopes A(t/R^2)^2 and |A|^2,
whose derivatives are O(R^-2), gives a contribution C_S R^-1 to both
c_R and c_R'. The error contributes at most C_S R^-3 T<=C_S R^-1
by the bounded oscillator propagator. Therefore

    |c_R(t)|+|c_R'(t)|<=C_S R^-1.                     (EF8)

This includes the zero-data switch-on transient, not just stationary
susceptibilities. Since g is exponentially localized,

    ||-12P_e(K P_R c_R g)||L2_6<=C_S R^-3/2.

LE then gives local response size C_S R^-3/2 sqrt(1+t). Its squared time
integral is bounded by C_S R^-3(T+T^2)<=C_S R, proving the shape part of
EF2. The shape mode has not been declared dispersive or omitted.

## 5. Direct cubic response and amplitude bookkeeping

The exact packet has L2 norm at most one for its real part and bounded
sup norm C_S R^-1/2. Thus

    ||P_e(P_R^3)||2<=C_S R^-1.

The same global-input CT estimate gives

    integral_0^T |q_D|^2dt<=C_S T^2 sqrt(1+T) R^-2
                                                   <=C_S R^3.

Using the square-of-three-summands inequality proves

    integral |q(Z_3)|^2<=C_(S,g4)(R^4+R+R^3).

Scaling the INPUT packet by epsilon scales Y_R and c_R by epsilon^2,
and each response in EF1 by epsilon^3. At R=epsilon^-1 the respective
squared trace budgets are O(epsilon^2), O(epsilon^5), O(epsilon^3).
This establishes EF3 without requiring cancellation between the channels.
The coefficients may have either sign; the estimate is not a claim that
the individual channels dissipate modal energy.

## 6. What has and has not been closed

The cubic continuous-response budget is closed for the single TP packet
comparison and its FULL second-order odd response. Spatial tails,
outgoing history, shape response, projection and switch-on are retained.
It is not a proof that the actual nonlinear solution stays close to this
third-order hierarchy on a critical interval. It is also not the separate
three-channel modal action-error theorem.

The next task is stability of this packet-response estimate under a
quantified packet remainder or superposition, including mixed interactions.
The actual generated continuum has a history and spectral distribution
that must be controlled, not selected to match TP's Gaussian. A sum of
single-packet bounds alone does not bound cubic cross terms without
summability/interaction estimates. Preserve that missing bridge explicitly.

ML's O(epsilon) whole-continuum bound for the original preparation remains
the admitted result in this chain. No incoming packet, fitted g4, new
source, clock, physical law or empirical FIN validation is introduced.
