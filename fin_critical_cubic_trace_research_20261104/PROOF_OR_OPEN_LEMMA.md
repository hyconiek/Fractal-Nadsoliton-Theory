# CT: sharp critical-time trace budget for the unweighted cubic forcing class

Date: 2026-11-04. Conditional analytical research; single-assistant audit.
Status: ACTUAL_CUBIC_COMPONENT_BOUND_WITH_SHARP_CLASS_OBSTRUCTION.
No actual critical damping law is claimed.

## 1. Correct frontier and unchanged variables

The continuation search found CD's three-channel response and VD's
threshold-defect virial, both already present in the source commit.
They supersede ST as the next-action frontier. Their existing spectral
densities, work form and continuum coercivity are reused, not rediscovered.
The task attacked here is the first obligation in VD's
NEXT_TRACE_AND_MODAL_GATE.md: control the actual compact mean from its
source history, including the unweighted cubic accumulation.

Use precisely VD's free Darboux phase (w_psi,pi_psi), obtained from QP by
continuous projection and the unitary Darboux map. It satisfies

    w_t=pi+f1, pi_t=w_xx-m^2 w+f2,
    m=sqrt(4+w_parameter), 2<=m<=sqrt(23/4),
    q(w)=<chi,w>, chi(x)=1_{[-1,1]}(x)/sqrt(1+x^2).

The w_psi field above is not the scalar parameter w_parameter. The QP
phase momentum is not replaced by the derivative of its position.
The actual even-continuum initial phase is zero. Split the EXACT supplied
forcing into F_loc+F_ge3 as in QP/VD, and define their zero-data linear
Duhamel responses in this free frame. This is an identity along the
actual solution, not autonomous evolution of either component:

    q_psi=q_loc+q_ge3.

For every finite time the paid input is
||F_ge3(t)||_(H1 x L2)<=C epsilon^3. It does not provide spatial moments
or temporal phase cancellation of this forcing. Global actual phase
energy is O(epsilon), and the physical initial state is unchanged.

## 2. A finite-time observation bound without weighted input

Let U_m(t) be the FREE Klein-Gordon phase group, and use
E_m(a,b)^2=||a_x||_2^2+m^2||a||_2^2+||b||_2^2.
For every T>=0,

    ||q(U_m(.) (a,b))||_L2(0,T)
                  <=C(1+T)^(1/4) E_m(a,b).              (CT1)

Here q acts on the position component only. C is uniform on the compact
mass interval. No positive spatial moment of (a,b) is assumed.

Proof for T>=1. Set delta=T^-1/2 and split Fourier space into
|k|<=delta and |k|>delta. With nu(k)=sqrt(k^2+m^2), the solution position is
cos(t nu) a_hat+sin(t nu)b_hat/nu. Since chi is in L1,
|chi_hat|<=||chi||_1. Cauchy-Schwarz on the low interval gives

    sup_t |q_low(t)|<=C delta^(1/2) E_m(a,b),
    ||q_low||_L2(0,T)^2<=C T delta E_m(a,b)^2.           (CT2)

For the high interval, split k into positive/negative half-lines and
the time factors into exp(+/-it nu). On each half-line nu is one-to-one.
After changing variables from k to nu, temporal Plancherel gives for
each component, up to fixed Fourier constants,

    ||q_high,branch||_L2(R_t)^2
      <=C integral_(|k|>delta) |chi_hat(k)|^2
                 (|a_hat(k)|^2+|b_hat(k)|^2/nu(k)^2)
                                      nu(k)/|k| dk
      <=C delta^-1 E_m(a,b)^2.                         (CT3)

The last bound uses 0<delta<=1 and the uniform mass band; nu/|k|
is at most C/delta. Sum the four components by Cauchy-Schwarz, then
restrict to (0,T). CT2 and CT3 give C sqrt(T) for the squared norm.
The sharp spectral cutoff is only a proof decomposition of L2 data;
no spatial H1 multiplier derivative is required. Density extends the
argument from Schwartz data. For T<=1 use global energy and chi in L2.
This proves CT1. In particular an all-time, unweighted-input uniform
observation constant has NOT been claimed.

## 3. Apply CT1 to the actual cubic component

Minkowski's inequality and time translation of CT1 give for zero-data
Duhamel with any integrable phase forcing F,

    ||q_Duh||_L2(0,T)
       <=C integral_0^T (1+T-s)^(1/4) E_m(F(s)) ds
       <=C(1+T)^(1/4) integral_0^T E_m(F(s)) ds.        (CT4)

An instantaneous phase impulse is (f1,f2); hence F's H1 x L2 norm is
exactly the regularity needed here. No derivative of f2 or of the source
history is taken. Apply the actual QP/VD forcing estimate:

    integral_0^T |q_ge3(t)|^2 dt
                           <=C epsilon^6 T^2 sqrt(1+T). (CT5)

For every fixed S<infinity and 0<epsilon<=1, this yields

    integral_0^(S epsilon^-2)|q_ge3|^2 dt<=C_S epsilon. (CT6)

This improves the O_S(1) bound obtained by squaring and integrating the
old pointwise C epsilon^3 t estimate. It is a bound on the cubic response
driven by the ACTUAL solution's forcing, but not on the entire trace.
The exact decomposition implies

    integral_0^T |q_psi|^2
       <=2 integral_0^T |q_loc|^2+C epsilon^6 T^2 sqrt(1+T). (CT7)

Inserting CT7 into VD's unchanged continuum virial gives on critical time

    integral E_rad^2 <=C_S epsilon^2
                              +C integral |q_loc|^2+C_S epsilon. (CT8)

This does NOT close the desired O(epsilon^2) trace budget. The purpose
of the next section is to decide whether the gap is merely a loose
estimate for this norm class. It is not.

## 4. Sharpness with zero initial data and an O(epsilon) energy bound

This is an ANALYSIS COUNTERCONTROL for the forcing norm class. It is NOT
the RP trajectory, a substituted QP forcing, an added FIN bath, or a new
physical initial packet. It tests an inference from bounds alone.

Fix any m in the same mass interval, 0<epsilon<=1/4, and set

    T=epsilon^-2, L=epsilon^-1, s=t/T, y=x/L,
    h(y)=1/(1+y^2), a(s)=s^2 for 0<=s<=1,
    W_epsilon(t,x)=epsilon L^-1/2 a(s)h(y)cos(mt),
    Pi_epsilon=partial_t W_epsilon.

Both initial coordinates are exactly zero because a(0)=a'(0)=0.
Choose f1=0 and DEFINE the comparison forcing by the exact residual:

    f2=epsilon L^-1/2[
       T^-2 a''(s) h(y)cos(mt)
       -2m T^-1 a'(s)h(y)sin(mt)
       -L^-2 a(s)h''(y)cos(mt)].                       (CT9)

Then W_tt-W_xx+m^2 W=f2 identically. The mass terms cancel, not merely
approximately. All functions are real, even in x and smooth for t>=0.
Since ||L^-1/2 h(x/L)||_2=||h||_2 and each x derivative costs L^-1,

    sup_(0<=t<=T) E_m(W_epsilon,Pi_epsilon)<=C epsilon,
    sup_(0<=t<=T)||f2(t)||_2<=C epsilon^3.               (CT10)

The first bound includes the physical velocity and spatial derivative;
the second uses T^-1=L^-2=epsilon^2 and T^-2=epsilon^4.
There is no localized quadratic forcing. One may smoothly extend a(s)
beyond s=1 to vanish for s>=2 with bounded derivatives; the same global
energy/source bounds then hold for all finite t. The stated interval
already suffices to test the claimed critical estimate.

The compact mean is EXACTLY

    q(W_epsilon(t))=epsilon^(3/2) A_epsilon s^2 cos(mt),
    A_epsilon=integral_-1^1 chi(x)h(epsilon x)dx.

For 0<epsilon<=1, 1/2<=A_epsilon<=2: on [-1,1], both
chi and h(epsilon x) are at least 1/2 and the interval has length 2.
Moreover A_epsilon -> 2arsinh(1). Therefore

    integral_0^T q(W_epsilon)^2 dt
       =epsilon A_epsilon^2 integral_0^1 s^4 cos^2(mTs)ds. (CT11)

Write cos^2=(1+cos(2mTs))/2. One integration by parts yields

    |integral_0^1 s^4 cos^2(mTs)ds-1/10|<=1/(2mT).

For mT>=10 this integral is at least 1/20. Our epsilon range gives
mT>=32. Thus, with explicit constants,

    epsilon/80 <= integral_0^T q(W_epsilon)^2dt
                                                <=4epsilon/5. (CT12)

For fixed m the ratio to epsilon tends to (2arsinh(1))^2/10>0;
convergence is also uniform in the displayed compact mass interval.
In contrast, the ratio to epsilon^2 diverges. This proves that O(epsilon)
in CT6 is the correct general order under the stated unweighted cubic
forcing assumptions, EVEN with zero data, correct even parity and a
uniform O(epsilon) phase-energy bound. It is not a lower bound for the
actual nonlinear source. The broad low-frequency profile is exactly
the structure that unweighted norms permit and positive moments exclude.

The unitary Darboux map transports this comparison phase into the even
continuous subspace of L with the same energy-order bounds. No new bound
mode is introduced. This strengthens its relevance as a class control;
it still does not make its forcing equal to the actual polynomial QP
functional or enforce the parent's conserved energy/action.

## 5. Precise failure and the remaining source-specific obligation

The attempted route was to obtain the required trace budget from the
existing actual-source norm bounds through a sharp free observation
estimate. CT6 is its positive result; CT12 is its adverse limitation.
No weight scan or forcing refit was used to manufacture closure.

The next missing input is a quantitative property of the ACTUAL transformed
cubic forcing that excludes or cancels the near-threshold accumulation
allowed by CT9. It must be derived from the QP polynomial source and
actual modal history. A zero threshold moment cannot be imposed: ET/CD
already show the even channel need not be threshold-null. Temporal
oscillations may provide the needed cancellation, but neither energy
smallness nor CR's frozen secular computation proves it.

The localized trace q_loc remains a separate unpaid history estimate.
The action-law error must also retain CD's THREE channels. CR's earlier
single-channel candidate is only a provisional diagnostic; it is not a
complete two-mode action closure. In particular the small initial shape
action must not be assumed to remain small on critical time.

Even if the whole continuum budget were only O(epsilon), an absolute
estimate for the leading mixed action error epsilon^2 integral E_rad
would give O(epsilon^(3/2)) by Cauchy-Schwarz on T=S epsilon^-2.
That is not o(epsilon^2). A better trace bound alone is therefore not
permission to declare averaged damping; the phase-resolved action error
still has to be controlled. No actual instability or universal no-go is
inferred from this stopped comparison-class result.

All supplied-model, ontology, strict/legacy, selector, physical units,
source and F01/F03/F06/ToE boundaries remain unchanged.
