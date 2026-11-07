# TP2 - uniform spatial kernels close the conditional cubic trace budget

Date 2026-11-07. Parent 49cdac4318017d8d61c7c73b1c91944a2bd3785f.
Self-audited conditional analysis of the supplied RP model.

## 1. Separate the threshold value of the complete actual source

Use CE's complete converted real even cubic phase forcing
F_3=(F_1,F_2) in the free Darboux frame. CE proves its global energy norm
is at most C epsilon^3 and its low-frequency energy density is at most
C epsilon^3. Its value at k=0 is gamma(t)=(gamma_1(t),gamma_2(t)) from CE8.
The L1 pairings in CE imply continuity in k for each fixed time by dominated
convergence. They also give |gamma(t)|<=C epsilon^3 throughout the critical
window. No time derivative of gamma is needed below.

Fix kappa in (0,1]. Choose a real even smooth cutoff chi supported in
|k|<kappa, equal to one on |k|<=kappa/2. Define phase sources by their
unitary Fourier transforms

    Fourier(F_ref)(t,k)=chi(k) gamma(t),
    F_rem=F_3-F_ref.                                  (TP2.1)

The reference spatial profile is a fixed Schwartz function. It is a
linear decomposition of the already converted actual forcing. The
zero-data responses satisfy the exact identity q_3=q_ref+q_rem. The two
phase components of gamma are kept separately.

## 2. The entire threshold history is controlled without a modulation bound

Let nu(k)=sqrt(m^2+k^2) and let q use the compact observation kernel of CT.
The scalar kernels for the two reference phase components are

    K_1(t)=integral conjugate(q_hat(k))chi(k) cos(t nu(k)) dk,
    K_2(t)=integral conjugate(q_hat(k))chi(k) sin(t nu(k))/nu(k) dk.

Both obey |K_l(t)|<=C(1+|t|)^-1/2, uniformly in the original mass band.
Here q_hat is smooth with bounded derivatives on the compact support,
since the spatial observation kernel is integrable and compactly supported.
For |t|>=1, nu''=m^2/(m^2+k^2)^(3/2) has a positive uniform lower bound
on that support. The one-dimensional oscillatory integral bound follows
by splitting at |k|=|t|^-1/2. The inner interval contributes O(|t|^-1/2).
On each outer half interval integrate exp(it nu) by parts using
|nu'|>=c|k|. The boundary at |k|=|t|^-1/2 and the integral of
|t|^-1/k^2 are O(|t|^-1/2). The amplitude derivative term is bounded
by C log(2+|t|)/|t|, which is also O(|t|^-1/2). Smooth cutoff endpoints
supply no additional term. Small time is bounded by the amplitude L1 norm.

The exact reference response is

    q_ref(t)=integral_0^t [K_1(t-s)gamma_1(s)
                                    +K_2(t-s)gamma_2(s)] ds.

Using the actual amplitude bound from CE gives

    |q_ref(t)|<=C epsilon^3 sqrt(1+t),
    integral_0^T |q_ref(t)|^2dt<=C epsilon^6(T+T^2).

Thus at T=S epsilon^-2,

    integral_0^T |q_ref(t)|^2dt<=C_S epsilon^2.        (TP2.2)

This includes resonant temporal histories, arbitrary measurable modulation
within the amplitude bound, and the source switch-on at zero. In particular
control of gamma's oscillation or its derivative is unnecessary for this
component estimate. This pays the leading threshold profile directly. Section 5 below uses
the stronger endpoint structure of CE to control the full residual as well.

## 3. A sufficient spectral-variation criterion

Write e_rem(t,k)^2=nu(k)^2|Fourier(F_rem,1)|^2
                                  +|Fourier(F_rem,2)|^2.
The global source energy and its low-frequency density remain bounded
by C epsilon^3. Apply TB's core/intermediate estimate and FH's high-output
estimate to F_rem. The core |k|<=epsilon already costs at most C_S epsilon^2
in squared trace because its energy density is bounded and its width is
epsilon. The transition region kappa/2<=|k| and all higher output also cost
O(epsilon^2) for fixed kappa. Therefore

    ||q_rem||L2(0,T)
       <=C_S epsilon
         +C integral_0^T [integral_(epsilon<|k|<kappa/2)
                                  e_rem(s,k)^2/|k| dk]^(1/2) ds. (TP2.3)

Together with CS4, a sufficient condition for the actual whole-continuum goal is the last
source-time integral being O(epsilon). The logarithm in CE's upper bound
can now be assigned to variation of the source around its threshold value.

A sufficient, stronger condition is the uniform squared Dini estimate

    sup_(0<=s<=T) integral_(0<|k|<kappa/2)
                                  e_rem(s,k)^2/|k| dk<=C epsilon^6. (TP2.4)

For example a uniform Holder bound

    |Fourier(F_3)(s,k)-gamma(s)|<=C epsilon^3 |k|^alpha,
    alpha>0,

implies TP2.4. The reference cutoff equals one in this region. Unlike the
zero-moment route discussed in TB, this sufficient condition permits a
nonzero gamma. A controlled first spatial moment of the complete converted
source is another sufficient condition, through the Fourier difference
estimate. The size bounds alone leave uniform moment or Holder control undetermined.
Section 5 supplies a direct route using the stronger CE preimage structure.
The time-integrated criterion TP2.3 also allows weaker hypotheses than
the displayed time-uniform sufficient condition.

## 4. Pointwise continuity and source size do not pay this modulus

A comparison family shows the distinction between continuity and a uniform
weighted difference estimate. Let g be a fixed even Schwartz function
with unitary Fourier transform g_hat(k)=exp(-k^2/2). Put

    G_L(x)=[g(x-L)+g(x+L)]/2,
    G_L_hat(k)=cos(kL)exp(-k^2/2), L>=1.

Its L1 and H1 norms are bounded independently of L, it is even, and its
threshold value equals one. In 0<k<kappa/2 its residual is
cos(kL)exp(-k^2/2)-1. The weighted square integral grows as

    integral_0^(kappa/2) |cos(kL)exp(-k^2/2)-1|^2 dk/k
                                       =(3/2)log L+O_kappa(1). (TP2.5)

For a proof first replace exp(-k^2/2) by one. The replacement changes the
integral by a bounded quantity since |exp(-k^2/2)-1|<=Ck^2. With y=kL,
(cos y-1)^2=3/2-2cos y+(1/2)cos(2y). Its integral divided by y is finite
near zero. On [1,kappa L/2] the cosine terms have bounded integrals by
one integration by parts, while the constant contributes (3/2)log L.
This proves TP2.5.

This is a forcing-class comparison. It shows why CE's size bounds and
pointwise continuity alone leave the remaining modulus unpaid. It supplies
no counterexample to the actual prepared nonlinear trajectory or to its
possible temporal cancellations. The first spatial moment of this family
grows with L, which is exactly information absent from the size bounds.

## 5. The stronger CE structure controls the whole residual directly

The complete physical phase source before projection, as supplied by CE3,
has the structure

    F_phys=(A_C, Omega_psi B_0+C_D),
    A_C=-(a/2)K(S_+-S_-),
    B_0=-(a/(2i))K(S_++S_-),
    C_D=-6v^2u-g4u^3,
    ||A_C||L1+||B_0||L1+||C_D||L1<=C epsilon^3.       (TP2.6)

CE1--CE2 pay the first two norms. TB6 pays the third. Its phase projection
and Darboux transform give the exact free F_3. Its global free phase norm
is already O(epsilon^3) by CS. Thus both an L1 preimage structure and a
global energy bound are available for the same complete source.

Let Q be the compact observation kernel, and retain the smooth fixed
cutoff chi from section 1. Set

    f_t=chi(D) cos(t Omega_psi) Q,
    g_t=chi(D) sin(t Omega_psi) Omega_psi^-1 Q.

For an impulse with phase source TP2.6 the low-output observation equals

    <Wcal f_t,A_C> + <Omega_psi Wcal g_t,B_0>
                                           +<Wcal g_t,C_D>. (TP2.7)

This follows by adjointness in the stated order. Physical projection
is accounted for by Wcal-adjoint P_e=Wcal-adjoint. In particular the
middle term keeps Omega_psi to the left of Wcal. It assumes no commutation
between these operators. The formula can first be read on smooth data
and extended using the bounds that follow.

For every integer j from zero to four,

    ||partial_x^j P^-1/2 f_t||Linfinity
       +||partial_x^j P^-1/2 g_t||Linfinity
                                       <=C_j(1+|t|)^-1/2. (TP2.8)

Here P=(-partial_x^2+1)(-partial_x^2+4) is the Darboux normalization.
Indeed each kernel is an integral of exp[i(kx +/- t nu(k))] times a
smooth compactly supported amplitude. The factor k^j, the normalization
multiplier, Q_hat and the optional nu^-1 are harmless on this fixed
support. The second phase derivative is +/-t nu''(k), whose magnitude
is at least c|t|. To see uniformity in x, split where the first derivative
has magnitude at most sqrt(|t|), an interval of length O(|t|^-1/2)
because it is monotone. On its complement integrate by parts. Boundary
terms and the variation of the inverse first derivative cost O(|t|^-1/2),
and the amplitude derivative has the same bound. This proves TP2.8
uniformly for all x, including |x| proportional to t. Small time follows
from the fixed amplitude L1 norm.

Write Wcal=T P^-1/2 with
T=partial_x^2-3K partial_x+3K^2-1. All coefficients of T and their first
two derivatives are bounded. Applying TP2.8 through order four gives

    ||Wcal f_t||Linfinity+||Wcal g_t||Linfinity
       +||partial_x^2 Wcal g_t||Linfinity
                                            <=C(1+|t|)^-1/2.

Finally use the exact identity

    Omega_psi Wcal g_t
          =Omega_psi^-1(-partial_x^2+m^2)Wcal g_t.

The Omega_psi^-1 convolution kernel has L1 norm 1/m by CE's heat formula.
It follows that all three test kernels in TP2.7 obey the same uniform
Linfinity bound C(1+|t|)^-1/2. Therefore the impulse observation satisfies

    |q_low(U_m(t)F_3(s))|
                    <=C(1+|t|)^-1/2 epsilon^3.       (TP2.9)

This controls the spatially extended source with its actual L1 preimages.
It requires neither a spatial moment nor the Dini condition TP2.4. The
translated comparison in section 4 illustrates why that sufficient Dini
criterion was stronger than needed for the available endpoint structure.

## 6. The conditional actual critical trace estimate

Integrating TP2.9 over the entire source history gives

    |q_3,low(t)|<=C epsilon^3 sqrt(1+t),
    ||q_3,low||L2(0,T)^2<=C epsilon^6(T+T^2).

For the complementary output multiplier 1-chi, FH9 applies because its
support is separated from zero by kappa/2. Its bounded multiplier size
and the actual global phase forcing norm give

    ||q_3,high||L2(0,T)^2<=C_kappa T^2 epsilon^6.

Add the two components by the square-of-two-summands inequality and set
T=S epsilon^-2. This proves, conditionally on the predecessor estimates,

    integral_0^(S epsilon^-2) |q_3(t)|^2dt<=C_S epsilon^2. (TP2.10)

CS4 then gives the actual whole local-continuum estimate and the whole
compact trace estimate

    integral_0^(S epsilon^-2) E_local(t)^2dt<=C_S epsilon^2. (TP2.11)

Constants depend on the fixed cutoff, S, finite model coefficients and
original mass band. The result includes CE's full cascade and its physical
phase conversion. It does not assume a zero threshold coefficient or a
special temporal carrier of the actual source. The paid inputs are the
actual global small-energy bounds, QP's bounded symbol derivatives and
CN/CS's exact decomposition.

## 7. Next research obligation and review status

The cubic trace gate is now conditionally closed by TP2.10--TP2.11.
The next task is the separate three-channel modal-action error estimate
using this critical local-continuum control. Reopen the exact CD/VD
obligation before choosing the error integrals to estimate.

These are self-audited analytical proofs in the supplied RP model.
Independent proof review and formal verification are outstanding.
A modal damping law still requires the action-error argument. The work
supplies no new physical preparation or empirical validation.
