# CR: source-bound secular coefficient and critical-time defect

Date: 2026-11-04. Base commit: 1dcba969e12fabbfd43695f800687f5c93076b75.
Status: PARTIAL_CONDITIONAL_REDUCTION; ACTUAL_CRITICAL_MODULATION_OPEN.
Single-assistant analytical work, not independent review or formal verification.

## 1. Inputs, scope and the one attempted route

Use the unchanged RP real odd/even Cauchy problem, w in [5/4,7/4],
omega=sqrt(w), ell=3, lambda=sqrt(ell), and finite signed g4,g6.
Reuse CM's exact canonical map and energy-space C1 remainder; SP's
global small-energy bound; ST's local propagation; and NM/OM's actual
memory kernel and causal primitive. No source, parameter or preparation
is changed. All claims below are conditional on those cited analytical
dependencies. No new independent verification of their PDE proofs is claimed.

The attempted route is to extract the resonant cubic feedback and compare
the actual canonical equations on T=S epsilon^-2 by a cumulative-defect
estimate. One refinement replaces a bound on the absolute forcing integral
by a bound on its oscillatory primitive. This still leaves an unpaid
estimate, explicitly identified in section 5; no second campaign is started.

New dependency delta: a cross-chart, normalization-explicit secular
coefficient for the frozen causal response, an exact actual-modal identity,
and a critical-time stability criterion with the fifth-order energy
remainder already paid. This is neither a new kernel nor a rediscovery of
NM's Volterra equation, OM's test-envelope theorem or CM's normal form.

## 2. Frozen causal feedback: a secular primitive, not a PDE envelope law

Let F=K f^2 and k(t)=<F,sin(t sqrt(H))/sqrt(H) F> on the odd sector,
INCLUDING the shape eigenfunction. In NM's normalization

    k=k_s+k_c, k_s=d^2 sin(lambda t)/lambda,
    chi(0)=1/70, Im chi(2omega)=Gamma/9,
    chi_s(nu)=d^2/(ell-nu^2), d^2=27pi^2/8192.

Gamma is the RP quantity for source 3F, not F. Write chi=chi_s+chi_c.
For a constant complex b define the algebraic diagnostic

    z_b(t)=(b exp(-i omega t)+conj(b)exp(i omega t))/2,
    N_b(t)=72 z_b(t) integral_0^t k(u)z_b(t-u)^2 du
                             -(18g4/35) z_b(t)^3,
    D_b(t)=(i/omega) exp(i omega t) N_b(t).

This evaluates the existing causal feedback functional on a frozen linear
carrier. It does NOT replace the actual nonlinear trajectory, and supplies
no external drive to the PDE. The original physical preparation corresponds
to real b=epsilon only in this diagnostic.

Then, uniformly over the frozen w interval and all t>=0,

    integral_0^t D_b(s) ds = t kappa |b|^2 b + O(|b|^3),       (CR1)
    kappa = -Gamma/omega
          + i[9(2chi(0)+Re chi(2omega))/omega
                                      -27g4/(140omega)].     (CR2)

The O constant depends on fixed g4 bounds and the inherited kernel, not t,b.

Proof. The three Fourier coefficients of z_b^2 at frequencies
0,-2omega,+2omega are |b|^2/2,b^2/4,conj(b)^2/4. For a constant envelope
OM's finite-interval identity is EXACT:

    convolution[e^(-i nu t)] = e^(-i nu t)[chi(nu)-P_nu(t)].

The stationary part of z_b times this response has coefficient
|b|^2 b [2chi(0)+chi(2omega)]/8 at -omega. The corresponding
coefficient of z_b^3 is 3|b|^2 b/8. Multiply by 72 and i/omega,
respectively, to get CR2. Realness implies chi(-nu)=conj(chi(nu)).

All other stationary demodulated terms have nonzero frequencies among
2omega,-2omega,4omega and therefore bounded primitives. The continuum
initial terms are absolutely integrable, since OM proves uniformly
P_c,nu in L1 for nu in {0,+2omega,-2omega}. The retained shape terms
e^(-i nu t)P_s,nu(t) have frequencies +/-lambda. After multiplication
by z_b and demodulation their frequencies belong to
{+/-lambda,2omega+/-lambda}; none is zero. Indeed
2omega-lambda >=sqrt(5)-sqrt(3)>1/2. Their primitives are bounded,
although the terms themselves need not decay. This proves CR1 without
discarding switch-on, shape or negative-frequency response. At t=0 the
exact convolution is zero, not its stationary susceptibility.

The coefficient is independent of g6 only at this cubic order; signed g6
remains in the exact dynamics and the higher energy remainder.

## 3. Cross-chart check and the continuum part that CM cannot remove

CM writes h=Z4=Az Jz^2+Bz Jz Jc+Cz Jc^2, where
Jz=P^2+w Z^2, Jc=R^2+ell C^2, and

    Az=27g4/(560w^2)+9d^2(8w-3ell)/(4w^2 ell(ell-4w)).

The pure-bound coefficient for a full-amplitude carrier b is
-4i omega^3 Az. Direct rational algebra gives

    9[2chi_s(0)+chi_s(2omega)]/omega-27g4/(140omega)
                                                     =-4omega^3 Az. (CR3)

Thus the frozen feedback's shape and self-cubic coherent shift agrees
exactly with CM. The part not encoded by that pure-bound energy is

    kappa_c=-Gamma/omega
             +9i[2chi_c(0)+Re chi_c(2omega)]/omega,
    kappa=-4i omega^3 Az+kappa_c.                         (CR4)

This is a normalization check across two descriptions of the supplied
model; it is NOT permission to identify CM and QP states. In particular
using Gamma for F rather than 3F, or treating b as half-amplitude, changes
the damping factor incorrectly. In half-amplitude A=b/2 the cubic damping
coefficient would instead be 4Gamma/omega.

CR1 shows that the retained resonant term has an order-epsilon primitive
on t=S epsilon^-2 when b=epsilon. It cannot be assigned to a uniformly
O(epsilon^2) error. This is a statement about the frozen diagnostic; no
counterexample to an actual full-PDE approximation is inferred from it.

## 4. Exact actual canonical identities and a paid critical remainder

Use the ACTUAL CM state (Z,P,C,R) and the physical continuous coordinates,
on which CM acts as identity. Let V3=E3m, V4=Q4m and V5=R5 in the EXACT
energy E2+h+V3+V4+V5. Define full-amplitude interaction variables

    b_z=e^(i omega t)(Z+iP/omega),
    b_c=e^(i lambda t)(C+iR/lambda),
    Tz(V)=e^(i omega t)(partial_P V-i partial_Z V/omega),
    Tc(V)=e^(i lambda t)(partial_R V-i partial_C V/lambda).

These use canonical momenta, not an assumed equality P=Z_t or R=C_t.
Hamilton's equations give EXACTLY

    b_z'=-2i omega(2Az Jz+Bz Jc)b_z + Tz(V3+V4+V5),
    b_c'=-2i lambda(Bz Jz+2Cz Jc)b_c + Tc(V3+V4+V5),      (CR5)
    Jz=w|b_z|^2, Jc=ell|b_c|^2.

For example, partial_P h=2P h_Jz and partial_Z h=2wZ h_Jz;
combining them gives -2i omega h_Jz b_z, fixing the sign.
The actual initial state is retained: b_z(0)=epsilon+O(epsilon^3),
b_c(0)=-d A_CM epsilon^2+O(epsilon^4), with zero canonical momenta.
At w=3/2 the shape slip vanishes exactly, but the physical forced shape
response remains. QP's different nonzero initial field shift is not reset.

The actual SP bound and CM inverse give state norm <=C epsilon for all
finite t. CM's C1 estimate therefore implies |T(V5)|<=C epsilon^4 and

    sup_(t<=S epsilon^-2)|integral_0^t T(V5) ds|
                                                  <=C S epsilon^2. (CR6)

This pays the fifth-order ENERGY remainder's effect at critical time in
modal integral equations. It does not pay V3 or V4. An energy O(epsilon^5)
bound alone would not suffice; the inherited derivative estimate is used.

If desired the coherent term can be removed by the EXACT phase
theta_z(t)=2omega integral_0^t(2Az Jz+Bz Jc)ds and its c analogue.
Then (e^(i theta_z)b_z)'=e^(i theta_z)Tz(V3+V4+V5).
These phases depend on the actual actions. They are O_S(1), not generally
o(1), on critical time. This identity does not prove slow C2 envelopes.

## 5. One cumulative-defect criterion and the first unpaid estimate

To test, rather than assume, a candidate critical reduction, define the
explicit cubic vector field G(b) by retaining both coherent terms in CR5
and adding (kappa_c |b_z|^2 b_z,0). This is a frozen-response candidate;
other shape/continuum channels may prevent it from being sufficient.
No actual averaged law is admitted here.

For the actual field solution define the exact defect

    D_epsilon(t)=integral_0^t [T(V3+V4)
                              -(kappa_c |b_z|^2 b_z,0)] ds. (CR7)

This definition includes all mixed channels, physical continuum positions
and momenta, actual initial transients, and the shape equation. There is
no residual declared negligible by naming it. In particular V3 contains
the actual source terms of ST, in physical continuum positions u,v:

    partial_Z V3=12 integral K(C g f u+Z f^2 v+f v u),
    partial_C V3=integral K(12C g^2 v+6g v^2+12Z g f u+6g u^2).

Let beta solve beta'=G(beta) from the exact actual CM initial coordinates.
Its coherent terms preserve both moduli; the real part of kappa_c is
negative, so |beta_z| is nonincreasing and |beta_c| is constant. Thus beta
stays O(epsilon) for all t>=0. On the common O(epsilon) ball, G has
Lipschitz constant C epsilon^2. Subtract its integral equation from CR5.
CR6 and Gronwall give, for fixed finite S,

    sup_(t<=S epsilon^-2)|b(t)-beta(t)|
       <=exp(C S)[sup_(t<=S epsilon^-2)|D_epsilon(t)|
                                                       +C S epsilon^2]. (CR8)

Only the cumulative signed/complex defect is used; no absolute-integral
smallness or derivative of its oscillatory forcing is required. Consequently
sup|D_epsilon|=o(epsilon) would suffice for critical modal approximation.
For full local field/physical-velocity approximation a separate whole-
continuum o(epsilon) estimate is also necessary; a modal comparison alone
does not supply it.

The FIRST UNPAID estimate in this route is precisely

    sup_(0<=t<=S epsilon^-2)|D_epsilon(t)|=o(epsilon).     (OPEN-CR)

What the inherited bounds actually yield: |T(V3)|<=C epsilon e(t),
|T(V4)|<=C epsilon^2 e(t), with ST
e(t)<=C[epsilon^2 sqrt(1+t)+epsilon^3 t+epsilon^2].
At critical time their absolute integrals are only O_S(1) and O_S(epsilon),
respectively. Neither proves OPEN-CR. The cubic candidate itself has an
O(epsilon) primitive. Its subtraction cannot be justified by triangle
inequality. CR1 computes the secular term for CONSTANT b, but applying it
pointwise to b(t) would discard actual envelope derivatives and other
mixed forcing. OM already exhibits why mere C1 test-envelope control is
insufficient for its shape estimate. QP free nonresonance does not remove
SP's nonzero mixed temporal resonances. The proposed cumulative criterion
is less demanding than a uniform C2 envelope theorem, but remains unpaid.

Therefore this route stops as PARTIAL, with neither an epsilon^-2 damping
law nor its negation proved. No claim of optimality for ST's 4/3 exponent,
all-time decay, physical light, natural action/geometry/clock, strict-kernel
membership, selector closure, or F01/F03/F06 closure follows.
