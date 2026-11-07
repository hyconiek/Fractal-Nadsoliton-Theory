# FH - causal packet interface and the actual spectral remainder

Date 2026-11-07. Parent 71d5f6d3b302e3a28d8b62f8b5974d5cb85bdcbb.
Self-audited conditional analysis of the supplied classical RP model.
The proof uses PS, PC, EF, CT and CS with their recorded hypotheses.

## 1. A causal representation in phase space

Let X_(a,theta) be the physical phase at time zero of a PS atom,
with common mass m, a in [a_-,a_+] and width aR. Write U_L for the
continuous even phase group and fix T=S R^2. Suppose the continuous
phase forcing has a strongly measurable representation

    F(s)=integral X_(a,theta) d mu_s(a,theta)+E(s),
    M_T=integral_0^T ||mu_s||TV ds < infinity,
    D_T=integral_0^T ||E(s)||energy ds < infinity.       (FH1)

The signed measures are real. Bochner integrability of the displayed
phase-valued integral is assumed. Phase space is H1 x L2 with the
positive continuous energy norm. Source position and momentum are
retained separately. For a physical momentum source, theta=pi/2 gives
X=(0,Wcal Omega_m h_(aR)), so (FH1) has concrete momentum-only examples.

The exact zero-data Duhamel response is

    V(t)=integral_0^t integral U_L(t-s)X_(a,theta)
                          d mu_s(a,theta) ds + V_E(t). (FH2)

Energy unitarity gives sup_t ||V_E(t)||energy<=D_T.
CT, conjugated by Wcal, gives

    ||q(Wcal-adjoint V_E,pos)||L2(0,T)
                           <=C(1+T)^1/4 D_T.           (FH3)

For T=S epsilon^-2, D_T=O(epsilon^3/2) is sufficient for the squared
trace of this linear remainder to be O(epsilon^2). Controlling its
nonlinear mixed interactions is an additional obligation.

## 2. Delayed packet births retain the cubic comparison estimate

For an atom born at s, use its phase U_L(t-s)X for t>=s and zero before
s. This is a phase impulse at s. The time derivative of its truncated
position may contain a delta term when the initial position is nonzero.
The equations below use the regular phase and ordinary time products of
positions. Their physical interpretation is the phase-forced system (FH1).

For each pair of delayed atoms, the quadratic source begins at
b=max(s_j,s_k). The odd response and shape oscillator start with zero
data at b. On [b,T] each input is a homogeneous solution of L, with
age in [0,T]. In the PC bilinear normal form, initialize the correction
at b and retain the term -U_H,c(t-b)Q[j,k](b). The estimates are uniform
in all birth times because every scaled age belongs to [0,S/a_-^2].
The EF dual local-energy argument then gives the same global odd phase
bound C and the shape amplitude bound C R^-1 as in PS4 and PS5.

For the shape calculation the two carriers are still +/-m. Time shifts
multiply their phases by constants of modulus one. Products have carriers
0,+/-2m, with C1 envelopes on [b,T]. Integration by parts retains the
boundary term at b, and 4m^2-3>=18 pays the uniform frequency separation.

A cubic outer atom can start before or after b. Its product with the
quadratic response is integrated only after all participating births.
The quadratic response keeps its entire evolution since b. PS7's source
bounds therefore give the same three trace norm bounds

    continuum feedback <=C R^2,
    shape feedback <=C R^1/2,
    direct cubic <=C R^3/2.                           (FH4)

Integrating the ordered triples against the source-time measures and
using Minkowski yields, for the packet part of FH2,

    ||Z_packet||Q<=C_(g4) M_T^3 R^2, R>=1.            (FH5)

First prove this for finite simple source-time measures. Their uniform
bounds and absolute integrability permit the Bochner limit. Finite-time
energy continuity of the quadratic and cubic Duhamel maps identifies
the limit with the full retarded triangular response. This keeps all
cross terms and every quadratic switch-on transient.

When M_T<=C epsilon and R=epsilon^-1, the squared cubic trace is bounded
by C epsilon^2. This condition concerns the total variation integrated
in source time. A bound on the energy norm of F gives a lower comparison
for its atomic norm, while an upper bound on that atomic norm requires
additional spectral structure. FH5 is conditional on that structure.

## 3. A concrete obstruction to a uniform whole-source packet fit

Set f=sqrt(3/4) sech^2(x) and phi=f^3. A physical even source component
is -g4 z(t)^3 P_e phi. For the prepared data, z(0)=epsilon. This example
concerns the nonzero coefficient g4 component separately. BT already
controls its localized response. It diagnoses the cost of fitting the
whole source before making that useful separation.

With the unitary Fourier transform and Darboux adjoint,

    Wcal-adjoint phi
      =P^-1/2 (partial_x^2+3K partial_x+2)phi,
    (partial_x^2+3K partial_x+2)phi=20phi-24sech^2(x)phi.

The Fourier recurrence

    Fourier(sech^8)=(k^2+36)/42 Fourier(sech^6)

follows by taking the transform of the second derivative of sech^6.
Consequently

    Fourier(Wcal-adjoint P_e phi)
        =-(4/7) sqrt((1+k^2)/(4+k^2)) Fourier(phi).    (FH6)

Projection drops out because Wcal-adjoint f=0. The unnormalized Fourier
integral of sech^6 equals

    pi k (k^2+4)(k^2+16) / [120 sinh(pi k/2)].

This formula follows from the sech^2 transform and two derivative
recurrences, and has its continuous value 16/15 at k=0. In particular
the phase G=(0,Wcal-adjoint P_e phi) has a strictly positive energy norm
c_I on the fixed annulus I={1/2<=|k|<=1}.

Each unit PS atom, at any age, has free spectral phase energy density

    (m^2+k^2) (aR/sqrt(pi)) exp[-(aRk)^2].

On I this gives the explicit bound

    ||1_I U_m(t)X_free||energy <= B_R,
    B_R=sqrt((m^2+1)a_+ R/sqrt(pi))
                                    exp[-a_-^2 R^2/8]. (FH7)

The total length of I is one. Thus a representation of G with atomic
mass M and annular residual norm d_I necessarily satisfies

    c_I<=M B_R+d_I.                                  (FH8)

For d_I<=c_I/2 the necessary mass grows at least as
c_I exp[a_-^2 R^2/8]/(2 sqrt((m^2+1)a_+R/sqrt(pi))).
This is a quantitative obstruction to uniformly bounded whole-profile
fits by these common-scale atoms. It applies to a normalized spatial
source component. Multiplying the source multiplies both sides by that
amplitude. The estimate alone makes no claim about accumulation of the
actual nonlinear history or the size of its final observation.

## 4. A paid actual component away from threshold

Use CS's exact cubic phase forcing after its physical projection and
Darboux conversion, denoted F_3 in the free even frame. Its paid estimate
is sup_t ||F_3(t)||energy<=C epsilon^3. For fixed 0<kappa<=1 define the
spectral projector P_hi onto |k|>=kappa on both phase coordinates.
It commutes with the free phase group. CT's half-line time Plancherel
calculation, now with this fixed cutoff, gives the all-time estimate

    ||q(U_m(.)P_hi X)||L2(R_time)
                         <=C kappa^-1/2 ||X||energy. (FH9)

Indeed the change from k to nu=sqrt(m^2+k^2) has Jacobian nu/|k|.
For |k|>=kappa it is at most C/kappa, and the Fourier transform of the
compact observation kernel is bounded. The two spatial half-lines and
two temporal branches are combined using Cauchy-Schwarz. High frequency
costs no additional derivative because the phase norm contains nu^2
on position and one power on squared momentum.

Apply Minkowski to the retarded integral and restrict each translated
full-line homogeneous estimate to [s,T]. This uses only integrable
phase forcing. Thus, on T=S epsilon^-2,

    ||q_3,hi||L2(0,T)<=C kappa^-1/2 T epsilon^3,
    integral_0^T |q_3,hi|^2dt<=C_S kappa^-1 epsilon^2. (FH10)

This is an actual cubic-component estimate, conditional on CS's paid
source bound. Constants are uniform in epsilon for each fixed kappa.
Choosing kappa proportional to epsilon restores an O(epsilon) bound,
so the intermediate shrinking bands still need attention.

With q_3=q_3,lo+q_3,hi, CS4 now implies

    integral_0^T E_local(t)^2dt
         <=C_(S,kappa) [epsilon^2+||q_3,lo||L2(0,T)^2]. (FH11)

The same reduction holds for the whole compact trace. The exact cubic
source includes the quadratic normal-form cascade terms of CS1.
The projector is applied after their full conversion to free phase.
It is applied to the output forcing, which avoids any assumption that
individual nonlinear input factors have separate spectral support.

## 5. Next research obligation and scope

The new dependency is FH10, which pays the actual cubic output at a fixed
distance from threshold. FH1--FH5 specify a causal packet route and its
source-time coefficient cost. FH6--FH8 explain why the full source needs
a residual or a spectral split before that route can be applied uniformly.

Next estimate the actual low-output-frequency forcing/history in dyadic
bands epsilon<|k|<kappa and the core |k|<=epsilon. Preserve its temporal
phases and the complete CS1 cascade. Seek a summable improvement over
FH9's band-edge cost, or establish a precise remaining obstruction.
Actual critical nonlinear closure and three-channel modal-action errors
are open. All statements concern the supplied model and the specified
comparison class. Numerical checks below test identities and quadrature.
Independent proof review and formal verification remain outstanding.
