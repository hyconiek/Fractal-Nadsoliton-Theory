# RP audit: bound mode and spectral density, not nonlinear radiation

Date: 2026-10-19. Conditional analytical mathematics, not mechanized proof,
independent peer review or physical validation. Incoming source is preserved.

## 1. Exact family and its restricted meaning

Interpret every displayed decimal parameter endpoint as its exact rational
value: Delta in[1/10,1/5], |d|<=1/100, |h|<=1/50, and a in[-3/4,-1/4].
For the stated family, the e1 channel decouples from the11-component
complement. Its scalar PT operator is -partial_x²+a+6-6sech²x. The normalized
mode sqrt(3/4)sech²x has eigenvalue a+2; the next channel level is a+5
and its essential threshold a+6. The complement form bound gives a strictly
higher ground level and a far-vacuum edge at least a+6.06. Thus the full
operator's essential edge is a+6 and its lowest real-block eigenvalue is
simple, with the two-real-dimensional complex phase doublet retained.

The source's strict bounds s_->399/200 and s_+<401/200 follow from rational
squaring; they pay complement separation>3/40 and the NS lower gap>49/40.
The source's decimal strings are numerical previews, not outward intervals.
These coarse rational margins suffice for the admitted theorem.

The family is **not open in the whole field-matrix space**: its many exact
zeros, including K_1j=0, define a constrained subfamily. Openness or robust
positive spectral density under arbitrary matrix perturbations is unproved.
It has the QW parameter *type*, not a selected legacy/strict kernel value.
Sparse K with a single2–3 edge is not the actual dense radial kernel.

## 2. The noncommuting block is a spectator here

Noncommutation occurs in channels2 and3, not in the chosen e1 mode. For
the displayed bound action, Psi self-interactions are component-radial,
Y is diagonal and K_1j=0. With every other component and velocity initially
zero, uniqueness keeps that entire complement zero in the full nonlinear
evolution. The spatial profile f need not remain fixed, but the **species**
subspace e1 is invariant. The active dynamics is therefore Phi plus one
complex Psi field, independently of h,d_2,d_3 and the other spectator data.

The exact Gamma is a useful conditional spectral witness. It is not a
kernel-sensitive signature of the noncommuting spectator block, nor a
theorem about arbitrary noncommuting active radiation channels. This
qualification does not invalidate the bound-mode or Gamma calculation.

## 3. Pay the scattering resolution, not only the Jost ODE

Let D_n=partial_x+n tanh x, D_n^*=-partial_x+n tanh x and H_f=-partial_x².
Then L_phi=D_2^*D_2, D_2D_2^*=H_1+3,
H_1=D_1^*D_1 and D_1D_1^*=H_f+1. With T=D_2^*D_1^*,

    L_phi T=T(H_f+4), T^*T=(H_f+1)(H_f+4).

Consequently T[(H_f+1)(H_f+4)]^-1/2, initially on a dense smooth core,
extends to an isometry from free L2 onto the continuous subspace of L_phi.
The orthogonal missing space is ker T^*, consisting exactly of the even
translation mode and odd shape mode: solving the two successive first-order
equations supplies those two L2 solutions and no others. This pays the
range statement. The source denominator (1-iq)(2-iq) differs from the
positive normalization sqrt((1+q²)(4+q²)) by a unit-modulus q-dependent
phase, which does not alter the spectral measure dq/(2pi).

Thus the source scattering states form the claimed continuum resolution.
Unit-modulus asymptotics or a few ODE residual tests alone would not have
paid completeness/normalization. No physical quantum interpretation is used.

## 4. Fourier recurrence, amplitude and both spectral roots

For F_n(q)=int sech^n x exp(-iqx)dx, differentiation/integration by parts gives

    F_(n+2)=(n²+q²)F_n/[n(n+1)],
    int tanh x sech^n x exp(-iqx)dx=-iqF_n/n.

For q>0, integrate sech²z exp(-iqz) around a lower rectangular contour of
height pi. The clockwise double pole at -i pi/2 has residue
iq exp(-pi q/2); the horizontal integrals differ by exp(-pi q). Vertical
parts vanish. Hence F_2=pi q/sinh(pi q/2), and the recurrence yields exactly
the source F_4 and F_6/F_4. No numerical quadrature establishes this identity.

Applying these identities to the conjugated Jost polynomial gives

    M(q)=(9/4) iq(q²+4)F_4/[8(1+iq)(2+iq)].

At energy Omega²=4+q² with q in[1,sqrt3], the spectral integral has two
simple roots. Its imaginary boundary value is

    (pi/(2pi)) sum_(r=+-q) |M(r)|²/(2q)=|M(q)|²/(2q).

This pays the factor of two and the retarded sign in Gamma. Exponential
source localization and the smooth nonthreshold Jost transform justify the
Lorentzian approximate-identity limit. Removing the two discrete modes
does not change this on-shell amplitude. The source Gamma formula is correct.

## 5. A coarse independent uniform positive certificate

No precision campaign is needed. On1<=q<=sqrt3, use3<pi<22/7,
sqrt3<7/4 and e<3. Then pi q/2<11/4<3 and
sinh(pi q/2)<27/2, so F_4>5/27. Also
81q(q²+4)/[2048(1+q²)]>=405/8192. Therefore

    Gamma > 125/73728 >0

uniformly on the frozen a interval. This is deliberately conservative and
depends on exactly the stated scattering/source normalization. Elementary
checks for pi>3 use the8-term integrated alternating geometric lower sum;
pi<22/7 follows from the positive integral x^4(1-x)^4/(1+x²); e<3 follows
from its factorial series and n!>=2^(n-1). These are not measured rates.

The incoming center value approximately0.106376851619 is a point diagnostic,
not an outward decimal certificate or emitted nonlinear power.

## 6. Fixed-time expansion and trace regularity

For the actual polynomial source family and smooth localized data, the
flow is smooth/analytic in epsilon on each fixed finite interval. Psi->-Psi
symmetry makes Psi odd and Phi even in epsilon. Differentiation of the
full PDE gives the displayed linear mode and second-order retarded scalar
equation, with zero independent scalar-wave data. Fixed-time Duhamel
estimates pay Psi=epsilon Psi_1+O_T(epsilon³) and
Phi=K+epsilon²u+O_T(epsilon⁴).

H1 x L2 estimates alone do not give pointwise derivative traces at +/-R.
For the displayed boundary-current expansion use the corresponding smooth
H2 x H1 (or higher) fixed-time estimates, which follow from the same local
polynomial wave argument for these smooth data, or use weak local-energy
balance. This upgrade still supplies **no uniform growing-time estimate**.

## 7. Preserve RP-004 as open

The order-epsilon² kink/bound-mode interference at fixed R can overwhelm
the putative order-epsilon⁴ radiative current. NS orbital stability and
positive Gamma do not pay its sign, initial slip, arrival or nonlinear error.
The source's exact stopped disposition ON_SHELL_COUPLING_ONLY is correct.
No nonlinear outward flux, relaxation, scattering or second-body absorption
is admitted. The next step must target that actual missing estimate, not
another amplitude decimal or noncommutation label.
