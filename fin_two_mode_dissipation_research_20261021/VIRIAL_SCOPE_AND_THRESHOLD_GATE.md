# Virial applicability, an actual obstruction, and its causal threshold repair

Date: 2026-10-21. Conditional same-author proof and primary-source scope audit.
No critical-time nonlinear virial theorem is asserted.

## 1. What the primary paper actually supplies

[Kowalczyk--Martel--Munoz, arXiv:1506.07420v2](https://arxiv.org/pdf/1506.07420v2), proves local asymptotic
stability for one scalar phi4 kink under ODD energy perturbations.
Its Theorem1.2 bounds the fourth power of its single internal mode and
weighted radiation energy integrated over all times. Sections2 and4--5
construct a virial/FGR functional coupling that mode to odd radiation.

The psi=0 invariant axis of the supplied RP model matches that scalar
equation after (t_KMM,x_KMM)=(sqrt2 t,sqrt2 x). Our ACTUAL preparation
has psi=epsilon f!=0 and includes another bound mode and EVEN radiation.
That theorem is not a theorem for this coupled trajectory.

The inspected proof also reports numerical comparisons/shooting values
in sections4.2--4.3. They are not a new outward FIN interval certificate
or a vector coercivity calculation. The present work uses the method as
a primary benchmark, not those numbers as payment of a different theorem.
The exact reading/provenance scope is in LITERATURE_REGISTER.json.

## 2. A concrete obstruction to the unchanged vector estimate

Let C=Pc,e(f^3), already an actual source profile in ET, and consider the
linear continuous phase solution of L=H+w with data(C,0).
Its threshold moment is <r0,C>=-4sqrt3/35!=0.
LE's exact rank expansion gives

    c log T-C0 <= integral_1^T ||U_m(t)(C,0)||X_-6^2 dt
                                        <=C1 log(2+T).   (1)

Indeed the leading2x2 oscillator matrix has no simultaneous zero in its
two phase coordinates and nonzero determinant m^2; its norm is bounded
below by c t^-1/2 for this moment vector. The remainder is O(t^-3/2).
The upper bound follows from LE's t^-1/2 estimate.
The lower bound already holds on a fixed bounded interval where r0 is
nonzero, so it also applies to a positive exponential local weight.

Thus an unmodified estimate of the form
integral_0^T ||radiation||local^2 <=C||initial||weighted^2,
with C independent of T and valid for all even continuous data, is false.
FGR positivity of the frozen two-mode action-work form does not cure this
linear threshold fact. A scalar odd-field virial cannot be imported
unchanged into the even channel.

This is a comparison-class obstruction, not a counterexample to the
special prepared nonlinear trajectory. It neither disproves critical
damping nor permits ignoring its threshold. ET separately proves the
actual cubic forcing is not identically spatial-threshold-null even
when g4=0. No unexplained orthogonality may be imposed as a repair.

## 3. One directly motivated repair: exact LOCAL threshold deflation

Fix a smooth analysis cutoff chi(t), zero for t<=1 and one for t>=2.
For the even continuum phase define the LE leading operator

    T_m(t) h = chi(t) r0(x)/sqrt(2pi m t)
            (m cos(theta)mu_a+sin(theta)mu_b;
             -m^2 sin(theta)mu_a+m cos(theta)mu_b),
    theta=mt+pi/4, h=(a,b), mu_h=<r0,h>.

Then E_m(t)=U_m(t)-T_m(t), with U_m already continuously projected as
in LE, obeys the derived operator bound

    ||E_m(t)||_(X_6->X_-6)<=C(1+t)^(-3/2).              (2)

This is exact subtraction of the known linear rank term in a LOCAL norm,
not addition/removal of physical radiation. r0 is not globally L2.
No global finite-energy change of variable is claimed by this subtraction.

QP supplies the ACTUAL localized forcing phase vector

    F_loc=( (R_loc^++R_loc^-)/2,
                  Omega(R_loc^+-R_loc^-)/(2i) ),

with ||F_loc(t)||X_6<=C epsilon^2. Define from the actual source history

    Theta(t)=integral_0^t T_m(t-s)F_loc(s)ds.             (3)

Theta has only two scalar causal coefficients multiplying r0. It is an
analysis functional, not a new microscopic state, detector or independent
source. Its coefficients depend on the already supplied PDE. It has the
a priori bound ||Theta(t)||X_-6<=C epsilon^2 sqrt(1+t), not a proved
small critical-time amplitude.

Using QP's exact phase Duhamel formula, its actual zero even-continuum
initial data, (2), and ||F_ge3||X_0<=C epsilon^3 yields

    ||Pc,e(psi,psi_t)(t)-Theta(t)||X_-6
                                  <=C(epsilon^2+epsilon^3 t). (4)

The physical inverse QP correction is O(epsilon^2) in global phase energy
and has been included. No derivative of a forcing history or slow-envelope
assumption enters this result. Equation(4) is a bound for the ACTUAL
prepared trajectory, though Theta is computed from its unknown history.

Thus all the nonintegrable LOCAL quadratic even response is isolated in
the explicit threshold history (3). The remaining global cubic response
can still be order epsilon at t=S epsilon^-2; it has not been paid away.

## 4. The first unpaid joint inequality

The critical proof needs a threshold-aware nonlinear estimate for
Theta together with the retained cubic/shape/modal terms. A possible
virial or modulated energy must include its actual source history and
retain all boundary terms. Its coercivity cannot be inferred from the
positive3-channel density vector alone.

The exact next proof target is a source-bound control of the time-integrated
error in replacing actual corrected action increments by the backreaction
work form in THREE_CHANNEL_RESPONSE. In particular it must pay both:

- the threshold-history contribution (3) on the full critical interval;
- the global cubic accumulation, which LE bounds only by C epsilon^3 t.

At t=S epsilon^-2 the current bounds permit order-epsilon local radiation
and order-epsilon-squared action changes. They do not show that the
error in an averaged action law is o(epsilon^2). That precise small
remainder, not more density digits, is still OPEN.

The scalar virial route was tested and its naive even extension failed.
Threshold deflation above is the one directly motivated repair executed
here, not a completed vector virial. Preserve the failure and this exact
scope. No guessed damping ODE, bath, clock refit or new physical law.

