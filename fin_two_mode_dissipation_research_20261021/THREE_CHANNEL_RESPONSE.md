# CD: all three quadratic radiation channels and their action-work form

Date: 2026-10-21. Conditional analytical REFERENCE-RESPONSE theorem.
This is not an averaged evolution theorem for the actual nonlinear PDE.

## 1. Exactly the supplied cubic interaction

Keep RP's unchanged real prepared sector, w=a0+2 in[5/4,7/4],
omega=sqrt(w), eta=sqrt3, H=-partial_x^2+4-6sech^2 x,
L=H+w, K=tanh x, f=sqrt(3/4)sech^2 x, g=sqrt(3/2)sech x tanh x.
In the actual cubic Hamiltonian the terms linear in radiation are

    E3,rad=6z^2<Fa,v_phi>+6c^2<Fb,v_phi>+12zc<G,v_psi>,
    Fa=K f^2, Fb=K g^2, G=K f g.                       (1)

Use their actual continuous projections for the reference waves. This
does not change any on-shell coefficient below. Discrete modes remain
separate. CM removes pure-bound cubic terms, not the interactions (1).
The coefficients6 and12 come from the supplied action, not symmetry alone.

Freeze trial amplitudes z=a cos(omega t), c=b cos(eta t) ONLY to define
the quadratic retarded-response coefficients. The reference forcing is

    v_phi,tt+H v_phi=-3a^2 Pc Fa[1+cos(2omega t)]
                         -3b^2 Pc Fb[1+cos(2eta t)],
    v_psi,tt+L v_psi=-6ab Pc G[cos((omega+eta)t)
                                      +cos((eta-omega)t)], (2)

with zero reference Cauchy data. These are analysis responses to terms
already in the source, not additional physical drives or initial waves.
A frozen trial does not solve the nonlinear parent indefinitely.

## 2. The three open channels and exact transforms

Use the already paid Jost convention D(k)=(1+ik)(2+ik), dk/(2pi),
S2=pi k/sinh(pi k/2), S3=pi(1+k^2)/[2cosh(pi k/2)],
S4=(4+k^2)S2/6. Direct Fourier recurrence gives

    M_Fa=3ik(k^2+4)S4/(32D),
    M_Fb=ik(k^2+4)(2-k^2)S2/(32D),
    M_G=[3/(2sqrt2)](1+k^2)(3+k^2)S3/(24D).              (3)

Fa and G are reused NM/ET profiles. The new Fb expression follows from
Fb=(3/2)[tanh sech^2-tanh sech^4] and
M_(tanh sech^p)=ik(p-1)(k^2+4)S_p/[p(p+2)D].

The open frequencies and momenta are

| channel | field | frequency | momentum squared |
|---|---|---|---|
| z^2 | phi |2omega|4w-4 in[1,3]|
| c^2 | phi |2eta|8|
| zc sum | psi |omega+eta|2sqrt(3w)-1 in(2,4)|

The zc difference frequency eta-omega is below the psi continuum:
(eta-omega)^2-(4+w)=-1-2sqrt(3w)<0. Both DC terms are below the odd
phi continuous spectrum. Their discrete responses are not resonant.
The two open phi frequencies are uniformly distinct; for example their
difference exceeds1/2. Thus no cross-frequency secular energy term
is hidden by combining them.

For a real localized profile J, operator A with far mass m and
nu^2=m^2+q^2, q>0, define the scalar density

    d_J(nu)=Im<J,(A-nu^2-i0)^-1 J>=|M_J(q)|^2/(2q).      (4)

The both-root normalization is essential. Equation(4) is a scalar
boundary value, not a bounded inverse on the continuous spectrum.

Denote the three densities by d_a(w),d_b,d_x(w). All are strictly
positive by (3). In particular this is a complete quadratic two-mode
radiative-channel test, not a finer numerical value of the old one.

## 3. Exact retarded mean work; the distinction from actual dynamics

For q_tt+Aq=-A0 J cos(nu t), zero data, the spectral response gives

    lim_(T->infinity) E_q(T)/T=(A0^2 nu/2)d_J(nu).        (5)

To pay the limit, the scalar response kernel and its modulated tail have
vanishing transients for every open frequency here (and the closed/DC
ones). Odd kernels decay as t^-3/2; even G has t^-1/2 decay but the fixed
carriers are separated from its threshold, so its oscillatory tail
vanishes. The same is true after the one time derivative needed for
work. The harmonic susceptibility therefore gives mean force-work
(A0^2 nu/2)Im chi. Different frequencies average to zero; the transient
contribution has zero Cesaro mean.

An independent factor check is the spectral oscillator identity: the
resonant part of cos(nu t) has coefficient1/2, energy has factor1/2,
and T^-1|integral_0^T exp(i(lambda-nu)s)ds|^2 tends to2pi delta(lambda-nu).
The two roots +/-q then give exactly (5), not twice or half that number.

Consequently the two reference waves in(2) have total mean energy growth

    P=9omega d_a a^4+9eta d_b b^4
                            +18(omega+eta)d_x a^2 b^2.   (6)

Closed/DC responses have zero linear-in-T energy growth. Frozen sources
would supply unbounded reference work; this is NOT energy produced by the
finite-energy parent at all times and not a newly proved emitted power.

Define I=omega a^2/2, J=eta b^2/2 as frozen canonical action LABELS, and
define the backreaction WORK functionals

    B_I=lim (1/T)integral_0^T (z_t/omega)F_z,rad dt,
    B_J=lim (1/T)integral_0^T (c_t/eta)F_c,rad dt,
    F_z,rad=-12z<Fa,v_phi>-12c<G,v_psi>,
    F_c,rad=-12c<Fb,v_phi>-12z<G,v_psi>.                 (7)

These B's are NOT derivatives of the constant trial labels I,J and NOT
yet actual dI/dt,dJ/dt. Direct trigonometric averaging gives

    B_I=-k_a I^2-k_x I J, B_J=-k_b J^2-k_x I J,
    k_a=36d_a/omega^2, k_b=36d_b/eta^2,
    k_x=72d_x/(omega eta).                              (8)

For example the mixed contribution to each action-work is
-18d_x a^2 b^2. The closed difference susceptibility is real and gives
zero average. No quantum-counting premise is used.

The work balance is exact in this reference calculation:
P=-omega B_I-eta B_J. It also follows by differentiating E3,rad:
reference wave work + bound trial work + interaction-energy derivative
is zero. The interaction term is bounded and has zero mean derivative.
Equation(8) supplies a candidate dissipative vector field ONLY after an
actual critical-time transfer proof, which is not provided here.

## 4. Coarse exact positivity bounds, not more Gamma digits

The old source has Gamma_old=9d_a and Gamma_old>125/73728.
Thus k_a=4Gamma_old/w>125/32256. The factor9 is retained.

At q_b^2=8, (3) gives |M_Fb(q_b)|^2=(3/8)S2(q_b)^2.
Since2<q_b<3,3<pi<22/7 and e<3,
pi q_b/2<5 and S2(q_b)>4/81. Therefore

    d_b>1/6561, k_b=12d_b>4/2187.                       (9)

For the mixed channel2<q_x^2<4 and S3(q_x)>1/18, since
cosh(pi q_x/2)<e^4<81. Equation(3) yields

    d_x=(1+q_x^2)(3+q_x^2)^2 S3(q_x)^2
                                     /[1024q_x(4+q_x^2)]
        >25/1769472,
    k_x>25/57344,                                      (10)

using omega eta<7/3. The elementary pi/e bounds are the same already
paid in the local RP spectral proof; no floating endpoint is promoted.

Since min(125/32256,4/2187)=4/2187 and I,J>=0,

    -(B_I+B_J)>= (2/2187)(I+J)^2.                       (11)

The mixed channel is positive as well, but is not needed to weaken
(11)'s simple uniform coercivity bound. These are exact reference-response
constants under the supplied couplings/clock, not natural constants or
proof of a coercive nonlinear virial functional.

## 5. Actual critical gate

In the actual prepared parent I is initially O(epsilon^2), while CM's
shape action is initially O(epsilon^4). Equation(8) by itself does not
prove those scalings persist or license autonomous action ODEs.
The full actual retarded history, localized threshold channel, cubic
terms and both initial slips still need the critical-time error bound.
The separate VIRIAL_SCOPE_AND_THRESHOLD_GATE explains the actual
obstruction to importing the scalar odd-field theorem unchanged.

