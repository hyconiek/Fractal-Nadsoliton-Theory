# BT: actual localized bound-source budget across an internal threshold

Date: 2026-11-04. Base commit: 0cbe925e7634b97690b76a922ec357ace492e87a.
Conditional analytical component theorem, not critical damping or a
new physical evolution law. Single-assistant audit; no formal verification.

## 1. Exact physical source split

Use the unchanged RP real odd/even solution and the original physical
projections, NOT the canonical CM positions:

    eta=c(t)g+v, psi=z(t)f+u,
    f=sqrt(3/4)sech^2 x, g=sqrt(3/2)sech x tanh x,
    H=-partial_x^2+4-6sech^2 x, L=H+w,
    w in[5/4,7/4], omega=sqrt(w), lambda=sqrt(3).

SP's actual small-energy estimate gives |c(t)|+|z(t)|<=C epsilon at
every finite time; it does NOT give c=O(epsilon^2) on critical time.
Let Pc denote the even continuous projection, removing f. The exact
physical equation is

    u_tt+Lu=Pc[-12K eta psi-6eta^2 psi-g4 psi^3-g6 psi^5],

with zero continuous Cauchy data. Define the bound-only cubic/quintic
source using the ACTUAL c,z:

    D(t)=-6c(t)^2 z(t) Pc(g^2 f)
          -g4 z(t)^3 Pc(f^3)-g6 z(t)^5 Pc(f^5).          (BT1)

Let u_D be its zero-data retarded response and u_R=u-u_D. The remainder
is the EXACT zero-data response to

    R=Pc[-12K eta psi
         -6(eta^2 psi-c^2 z g^2 f)
         -g4(psi^3-z^3 f^3)-g6(psi^5-z^5 f^5)].        (BT2)

No component of R is declared higher order or small without proof.
In particular it contains the bound-only QUADRATIC term -12cz Pc(Kgf),
as well as every continuum-bearing nonlinear term. This is an analysis
decomposition of the actual equation, not a finite-mode truncation.

## 2. A paid actual component estimate on critical time

The three profiles in BT1 are fixed, smooth and exponentially localized;
their continuous projections preserve every required polynomial weight.
For fixed finite bounds on signed g4,g6 and epsilon<=1,

    ||D(t)||_L2_6 <= C epsilon^3.                       (BT3)

Here L2_6 means ||<x>^6 D||_2. The strong version is also available
because spatial derivatives act on fixed profiles, not on unknown
remainder fields. No derivatives of c or z are required.

Use LE's proved phase propagator bound X_6->X_-6, with position in H1
and momentum in L2, on the continuous even sector. Apply its first-order
Duhamel formula to the phase forcing (0,D):

    ||(u_D,u_D,t)(t)||_X_-6
       <=C epsilon^3 integral_0^t (1+t-s)^(-1/2)ds
       <=2C epsilon^3[sqrt(1+t)-1].                   (BT4)

The velocity is the physical derivative here. Differentiating the sine
propagator introduces no endpoint term because sin(0)=0, and is already
covered by the phase estimate. BT4 is valid for arbitrary actual time
dependence of c,z subject to the paid amplitude bound. It assumes no
carrier separation, slow envelope or independent incoming wave.

Consequently for fixed S and T=S epsilon^-2,

    sup_(t<=T)||(u_D,u_D,t)||_X_-6 <= C_S epsilon^2,
    integral_0^T ||(u_D,u_D,t)||_X_-6^2 dt
                                             <= C_S epsilon^2. (BT5)

The pointwise bound follows from epsilon^3 sqrt(T)=O(epsilon^2).
For the integral, [sqrt(1+t)-1]^2<=t gives an upper bound
C epsilon^6 T^2. The quintic piece alone improves these to
C_S epsilon^4 pointwise and C_S epsilon^6 after squared integration.
The z^3 term reuses SP's existing U_4,z mechanism. The added payment is
the retained c^2 z term, quintic component, and complete trace accounting;
SP's earlier z^3 estimate is not claimed as a new discovery.

This is exactly the type of spatial structure missing from CT's generic
unweighted cubic forcing class. CT's broad comparison profile cannot be
identified with a fixed linear combination of these localized profiles.

## 3. Transfer to the actual VD compact mean, without identifying charts

Let Wcal* be VD's normalized Darboux map and q(a)=integral_-1^1
<x>^-1 a(x)dx. Define q_D=q(Wcal* u_D), q_R=q(Wcal* u_R) in the
physical chart. Weighted Darboux boundedness and compact support of q give

    integral_0^T |q_D|^2 dt<=C_S epsilon^2.             (BT6)

VD's q_psi instead observes the QP position. Reconstruct the physical
state using QP's exact inverse: its position correction has global H1
norm at most C epsilon^2 on the actual small-energy trajectory. Thus

    q_psi=q_D+q_R+delta_q,
    |delta_q(t)|<=C epsilon^2,
    integral_0^T |delta_q|^2 dt<=C_S epsilon^2.          (BT7)

This uses only position boundedness; QP's phase momentum is not replaced
by its position derivative. Its actual nonzero odd initial chart shift
and zero even-continuum initial data remain unchanged.

It follows that

    integral_0^T |q_psi|^2 dt
                      <=3 integral_0^T |q_R|^2 dt+C_S epsilon^2. (BT8)

VD's actual continuum budget therefore reduces to

    integral_0^T E_rad^2 dt
                               <=C_S epsilon^2+C integral_0^T |q_R|^2dt. (BT9)

This is a genuine removal of a controlled SOURCE COMPONENT from the
unpaid trace, not subtraction of physical radiation or a non-L2 state.
It does not assert a bound on q_R. In particular the quadratic cz source
and continuum-bearing products remain there, as explicitly listed in BT2.

## 4. Why a uniform nonresonance shortcut fails inside the original band

The BT4 route did not assume frequency separation. That matters because
the c^2 z product has a mixed carrier

    nu_*(w)=2lambda-omega,
    nu_*(w)^2-(4+w)=8-4sqrt(3w).                       (BT10)

It reaches the even continuum threshold m=sqrt(4+w) exactly at
w_*=4/3, strictly inside [5/4,7/4]. There m_*=4/sqrt(3).
This is not ST's sufficient time exponent 4/3; it is a parameter value.
The squared detuning has derivative -2sqrt(3)/sqrt(w), hence derivative
-3 at w_*. No uniform positive gap for this carrier exists on the band.

The direct source is not spatially threshold-null either. With
r0=(3tanh^2 x-1)/2, rational beta integrals give

    <r0,g^2 f>=sqrt(3)/35, <r0,f>=0.                  (BT11)

Indeed substitution y=tanh x leaves (3sqrt(3)/8) times
integral_-1^1 (4y^4-y^2-3y^6)dy=(3sqrt(3)/8)(8/105).
Continuous projection therefore does not remove BT11.

BT10-BT11 concern a term in the actual polynomial source. They do not
prove that its full nonlinear history is a monochromatic threshold drive.
The remaining source terms could in principle interfere. The following
frozen two-mode calculation checks the quadratic cascades as well, while
keeping that narrower scope explicit.

## 5. A nonzero dressed frozen threshold coefficient

This section is a REFERENCE-RESPONSE coefficient calculation, not an
actual critical-time expansion. Introduce independent full complex trial
amplitudes b,c for physical first-order modes:

    psi_1=(b e^(-i omega t)+conj(b)e^(i omega t))f/2,
    eta_1=(c e^(-i lambda t)+conj(c)e^(i lambda t))g/2.

Write Fb=K g^2, Gx=K g f, J=r0 K f. Define the full odd outgoing
response R_H=(H-4lambda^2-i0)^-1 and the full even sub-spectrum
response R_L=(L-(lambda-omega)^2)^-1. The first is a weighted local
boundary value, never a bounded inverse on global L2. The second is a
real bounded inverse: w-(lambda-omega)^2=2sqrt(3w)-3>0 on the band.

The second-order eta coefficient at 2lambda is -(3/2)c^2 R_H Fb.
The second-order psi coefficient at lambda-omega is -3c conj(b)R_L Gx.
Both come from the SAME source's quadratic terms, with causal harmonic
coefficients in the retarded response. Discrete switch-on oscillations
and continuum transients are retained in the finite-time response; they
are not part of its selected stationary harmonic coefficient.

Substitution into -12K eta psi-6eta^2 psi gives the coefficient of
c^2 conj(b)e^(-i(2lambda-omega)t) as the localized even profile

    A_*=9K f R_H Fb+18K g R_L Gx-(3/4)g^2 f.           (BT12)

The g4 psi^3 term has b^3 harmonics, not c^2 conj(b), and g6 starts at
degree five in these independent trial amplitudes. No term in the
quadratic-cascade calculation for this particular cubic monomial is
discarded. The pure-bound intermediate components are included by the
full resolvents; one must not silently replace them with Pc alone.

The latter two terms in <r0,A_*> are real. From the inherited ET/CD
transforms, with D(k)=(1+ik)(2+ik), S2=pi k/sinh(pi k/2),

    M_J=-i sqrt(3) k^3(k^2+4)S2/(64D),
    M_Fb=ik(k^2+4)(2-k^2)S2/(32D).

At q_s^2=8, their ratio is M_J/M_Fb=2sqrt(3)/3, real and positive
at both roots. In the paid dk/(2pi) resolution the cross-density is
therefore

    Im<J,R_H Fb>=(2sqrt(3)/3)d_b,
    d_b=|M_Fb(q_s)|^2/(2q_s)>0.

Consequently

    Im<r0,A_*>=6sqrt(3)d_b>0.                          (BT13)

This rules out cancellation of the threshold coefficient within the
complete stationary quadratic-cascade plus direct cubic calculation for
that monomial. It uses the existing CD d_b, not a new Gamma precision
calculation. It is independent of g4 and g6 at this cubic jet order.
It does NOT exclude cancellations involving other orders or actual
history-dependent remainders of the prepared trajectory.

The independent two-mode trial is an analysis device, not a new Cauchy
preparation. In the ACTUAL prepared solution the physical shape starts
at zero. CM's free shape slip is only order epsilon^2; at w=4/3 its
coefficient is -2d epsilon^2/7+O(epsilon^4). It is not licensed as an
order-epsilon freely prepared shape. Establishing that small hierarchy
uniformly through critical time remains part of the open problem.

## 6. Threshold response phase, bounded component, and the remaining gate

LE gives for any fixed localized complex profile A and its even Pc part

    S_L(t)Pc A = r0<r0,A> sin(mt+pi/4)/sqrt(2pi m t)
                                                 +O_H1_-6(t^-3/2).

For the zero-data REFERENCE forced response at precisely m,
convolution with exp(-imt) then yields

    U_A(t)= e^(i pi/4) r0<r0,A> exp(-imt)
                                     sqrt(t)/sqrt(2pi m)+O_H1_-6(1). (BT14)

The physical velocity has the same leading expression multiplied by -im
and an O_L2_-6(1) remainder. It follows by applying LE's cosine-propagator
expansion before time integration, not by differentiating an unspecified
error term. The derivative of sqrt(t) only contributes a bounded tail.

Proof: multiplying the leading sine by exp(imt) leaves the constant
coefficient i exp(-i pi/4)/2=e^(i pi/4)/2. Its t^-1/2 integral is
2sqrt(t). The 2m oscillation has a bounded primitive and the remainder
is integrable; short times contribute only a bounded term. This is a
complex amplitude convention, not a claim about actual nonlinear growth.

This also records a narrow phase correction to ET equation (8): under
ET's own displayed sine(mt+pi/4) convention, its leading factor should
be exp(i pi/4), not i. The order sqrt(t), non-L1 conclusions and scope
remain unchanged. Historical ET files/hashes are preserved; the correction
is documented here rather than silently rewriting old proof bytes.

An internal threshold resonance does not by itself make BT5 fail:
the absolute t^-1/2 convolution already allows sqrt(t) growth. A fixed
localized source of size epsilon^3 still gives local size O(epsilon^2)
and integrated squared size O(epsilon^2) on critical time. Thus the
failed uniform-nonresonance shortcut is repaired by an estimate that
uses actual amplitudes and proven spatial localization instead.

The remaining concrete obligation is q_R in BT2/BT8, together with CD's
full three-channel modal action error. Prove source-specific control of
the retained quadratic cz history and continuum-bearing nonlinearities.
Do not identify physical and QP forcing coefficients, assume the small
shape hierarchy, set a threshold moment to zero, or promote BT13/BT14
to an actual epsilon^-2 damping or instability theorem.

F01, natural source/geometry/clock, strict/legacy completion, selector,
physical photons and ToE boundaries remain unchanged.
