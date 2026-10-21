# ET: even continuum thresholds and actual mixed-coupling cancellations

Date: 2026-10-21. Conditional same-author analytical proof.
Not a full-PDE epsilon^-2 theorem, independent review or physical source.

## 1. Exact source interaction; no two-profile truncation

Use unchanged F_RP_20261019_PREDECLARED, a0 in[-3/4,-1/4],
w=a0+2 in[5/4,7/4], finite signed g4,g6, and its actual real/odd/even
preparation. Put K=tanh x, H=-partial_x^2+4-6sech^2 x, L=H+w,
f=sqrt(3/4)sech^2 x, g=sqrt(3/2)sech x tanh x. Write the OLD physical
projections eta=c g+v, psi=z f+u, v orthogonal g, u orthogonal f.
CM leaves continuum positions/momenta unchanged but reconstructs z,c from
the new coordinates. Its actual initial slip must remain in applications.

Let Pc,e remove f from the even sector, and d=<K f^2,g>. Define

    B=Pc,e(K g f)=K g f-d f,
    C=Pc,e(f^3)=f^3-(18/35)f.                            (1)

The exact even continuum equation, with both initial u data zero, is

    u_tt+L u=-12 c z B-g4 z^3 C+Pc,e R,
    R=-12K(c g u+z f v+v u)-6eta^2 psi
      -g4(3z^2 f^2 u+3z f u^2+u^3)-g6 psi^5.            (2)

The z f v term generally has the SAME cubic order as c z. R is not
automatically a higher-order error. This task retains it explicitly.

## 2. Exact transforms and the even threshold state

Use the paid Darboux/Jost resolution with dk/(2pi) and
D(k)=(1+ik)(2+ik). If S_p is the Fourier transform of sech^p, then

    S_(p+2)=(p^2+k^2)S_p/[p(p+1)],
    M_(sech^p)=(2-p)(1+k^2)S_p/[(p+1)D],
    M_(tanh sech^p)=ik(p-1)(k^2+4)S_p/[p(p+2)D].          (3)

Expand the conjugated Jost numerator
3tanh^2-1-k^2+3ik tanh and use integration by parts
integral tanh sech^p exp(-ikx)=-ik S_p/p to prove (3).
In particular M_(sech^2)=0, paying the continuous removal of f.

Besides the inherited S2=pi k/sinh(pi k/2),

    S1=pi/cosh(pi k/2), S3=(1+k^2)S1/2,
    S6=pi k(k^2+4)(k^2+16)/[120sinh(pi k/2)].

For S1, the lower height-pi rectangle has its simple pole at -i pi/2.
Its clockwise residue is2pi exp(-pi k/2); the opposite sech sign on the
lower segment gives (1+exp(-pi k))S1=2pi exp(-pi k/2). Extend by evenness
and the removable value at0. The recurrence proves S3,S6. Consequently

    M_B=[3/(2sqrt2)](1+k^2)(3+k^2)S3/[24D],
    M_C=-[3sqrt3/14](1+k^2)S6/D,
    b0=M_B(0)=3pi/(64sqrt2), c0=M_C(0)=-4sqrt3/35.        (4)

The bounded EVEN threshold state is
r0=(3tanh^2-1)/2=1-(3/2)sech^2. It is not L2. Since integral r0 f=0,
the projections in (1) do not change (4). This resonance is a property
of the supplied operator, not an additional state, field or noise law.

## 3. A rank-one t^-1/2 response matrix

Set m=sqrt(4+w), G1=B,G2=C, and

    K_ij(t)=<Gi,sin(t sqrt(L))/sqrt(L) Gj>.

Its continuum integral has density conjugate(M_i)M_j against dk/(2pi).
The numerators in (4) are real even and their denominator is shared, so
the cross density is real even. Uniformly in the original parameter band,

    K_ij(t)=v_i v_j sin(mt+pi/4)/sqrt(2pi m t)
                                      +O(t^-3/2), v=(b0,c0). (5)

Proof, including remainder: set k=y sqrt(2m+y^2), turning the frequency
into m+y^2. The transformed amplitude is even Schwartz with uniform
seminorms. Gaussian-regularized Fourier inversion and Taylor expansion of
exp(-i xi^2/(4t)) give its value-at-zero leading term. The error uses a
finite weighted Fourier moment, exactly as in NM's odd-profile argument.
An extra expansion gives the second oscillatory term with O(t^-5/2)
remainder. Model-time derivatives are proved by multiplying the spectral
amplitude by its smooth frequency, not by differentiating an unspecified
O term. Small times are bounded directly by the spectral formula.

Both diagonal kernels fail to be L1. The leading matrix is rank-one and
positive semidefinite, with negative off-diagonal coefficient. It is not
two independent noises. The old odd F-F kernel estimate cannot be copied
to this even channel.

## 4. Exact threshold-null direction and its stronger response bound

Define only an ANALYSIS basis change of the existing profiles:

    Cstar=C+kappa B, kappa=256sqrt6/(105pi).              (6)

Then M_Cstar(0)=0. Its even numerator vanishes quadratically, and

    M_Cstar(k)=sqrt3(49-4pi^2)k^2/420+O(k^3).              (7)

For an explicit check use unnormalized B'=sech^3-sech^5, C'=sech^6.
The null profile is C'+1024 B'/(105pi). Taylor expansion gives

    S6=(16/15)[1+(5/16-pi^2/24)k^2+O(k^4)],
    S3=(pi/2)[1+(1-pi^2/8)k^2+O(k^4)].

The null transform's k^2 coefficient is2(49-4pi^2)/315. Multiply by
3sqrt3/8 to obtain (7). It is nonzero since pi<22/7 gives
49-4pi^2>465/49>0. Thus the self-kernel Kstar is O(t^-5/2), with finite
absolute first moment. Its cross kernels with B,C are O(t^-3/2).
This does NOT license replacing the physical forcing by Cstar.

For the general matrix, away from |nu|=m define successively

    P_nu(t)=integral_t^infinity exp(i nu u)K_ij(u)du,
    Q_nu(t)=integral_t^infinity P_nu(u)du, chi_nu=P_nu(0).

These are ITERATED oscillatory tails; do not interchange them with an
ordinary integral (u-t)K(u)du. The two-term expansion after (5) and
integration by parts show P,Q bounded (indeed O((1+t)^-1/2)) uniformly
for bounded nu with |nu-m|,|nu+m|>=delta>0. Oscillatory leading terms
have nonzero frequencies nu+/-m; remaining tails are integrable.
P'=-exp(i nu t)K, Q'=-P. chi is a scalar/matrix Abel resolvent boundary
value, not a bounded L2 inverse on the continuous spectrum.

For A in C2([0,S]), 0<epsilon<=1 and t<=S epsilon^-2, put

    V(t)=integral_0^t exp(i nu u)K_ij(u)A(epsilon^2(t-u))du.

Two integrations by parts give EXACTLY

    V=chi A(epsilon^2t)-P(t)A(0)+R,
    R=epsilon^2[Q(t)A'(0)-Q(0)A'(epsilon^2t)]
        +epsilon^4 integral_0^t Q(u)A''(epsilon^2(t-u))du,
    R_t=-epsilon^2 P(t)A'(0)
        -epsilon^4 integral_0^t P(u)A''(epsilon^2(t-u))du.

Hence sup(|R|+|R_t|)<=C epsilon^2[||A'||+(1+S)||A''||].
Derivatives are in ORIGINAL model time. No C3 norm or vanishing causal
slip is assumed. At nu=m, A=1 instead gives for diagonal entries

    V(t)=i v_i^2 sqrt(t)/sqrt(2pi m)+O(1).                (8)

So the general matrix theorem cannot include that temporal resonance.

For the null self-kernel Kstar, A in C1 suffices for value AND model-time
derivative error O(epsilon^2||A'||), uniformly in EVERY real nu, including
+/-m. Indeed ||P||1<=integral u|Kstar(u)|du gives the value remainder.
For its derivative use
R_t=epsilon^2[V_(A')(t)-chi A'(epsilon^2t)], bounded by
2epsilon^2||Kstar||1||A'||, without differentiating A'.

The leading bound-shape c z and direct z^3 terms in (2) use carriers
+/-sqrt(w), +/-3sqrt(w), +/-[sqrt(w)+sqrt3], +/-[sqrt(w)-sqrt3].
Their distance from +/-m is >1/4. The rational bounds
11/10<sqrt(w)<4/3,17/10<sqrt3<7/4,9/4<m<12/5 give smallest coarse gap2/5.
This is NOT the full forcing spectrum: z f v contains continuous transient
frequencies and remains in the exact remainder (2).

## 5. Stronger cancellation in the ACTUAL mixed threshold moment

The full odd-field contribution to the even forcing's threshold moment
uses the fixed ODD profile J=r0 K f. Equation(3) gives

    M_J=-[sqrt3/2] i k^3(k^2+4)S2/[32D],
    M_J/(ik^3)->-sqrt3/16.                               (9)

Indeed J=(sqrt3/2)[tanh sech^2-(3/2)tanh sech^4], and the two odd
transforms combine into the extra k^2 factor. With F=K f^2, NM gives
M_F/(ik)->1/4. Therefore the cross density is
conjugate(M_J)M_F=-sqrt3 k^4/64+O(k^6). Its odd CONTINUUM kernel obeys

    K_JF,c(t)=-3sqrt3 sin(2t+5pi/4)/(32sqrt(pi)t^(5/2))
                                             +O(t^-7/2). (10)

The proof is the same Gaussian expansion, now with the first two even
amplitude coefficients zero. Its absolute first moment is finite.
The BOUND shape piece is not deleted: <J,g>=d/sqrt3, <F,g>=d, hence
K_JF,b=(d^2/3)sin(sqrt3 t). CM/OM must retain that undamped component.

This useful cancellation is not full threshold decoupling. Let eta2 be
the actual second Taylor coefficient, defined at all reference times by

    eta2_tt+H eta2=-6F cos^2(sqrt(w)t), zero Cauchy data.

The FULL cubic source (not its shape-only approximation) has threshold
amplitude

    h(t)=-12cos(sqrt(w)t)<J,eta2(t)>
                         -g4 cos^3(sqrt(w)t)c0.          (11)

Rational beta integrals give

    <J,F>=0, <J,HF>=-32sqrt3/385.

For a direct check HF=-9tanh(x)sech^4(x)+18tanh(x)sech^6(x).
Substitution t=tanh(x) turns both pairings into polynomial integrals;
no numerical stationary-phase estimate is used for these identities.

If g4!=0 then h(0)=4sqrt3 g4/35!=0. If g4=0, use
eta2''(0)=-6F and eta2''''(0)=6HF+12wF to get

    h(0)=h''(0)=0, h''''(0)=2304sqrt3/385>0.              (12)

Thus for EVERY fixed g4 the actual cubic source is not identically
threshold-null. Finite signed g6 enters at epsilon^5 and does not change
this cubic conclusion. The early cancellation <J,F>=0 also shows why a
shape-only replacement is wrong: it loses a real cancellation with the
continuum. Nonzero SPATIAL threshold overlap does not mean temporal
forcing at frequency m, secular growth or a proved large-time field tail.
Taylor coefficients agree with the actual PDE on fixed finite windows;
no uniform nonlinear epsilon^-2 expansion is asserted.

## 6. Exact remaining full-PDE gate

Paid: actual even response matrix, rank-one threshold tail, exact null
direction, causal modulation bound, extra mixed J-F cancellation, and the
actual cubic nondecoupling test. No parameters, data or action were fitted.

Still open is the spatial two-propagator interaction, schematically

    S_L(t-s) Pc,e[K f S_H(s-r) Pc,o F],

with its true modulated coefficients and all remaining nonlinear terms
on the growing interval. A scalar threshold moment, a two-profile matrix
or CM's finite energy jet does not pay that spatial estimate or actual
slow-envelope regularity. No averaged damping, natural clock, stochastic
law or source selection follows. NM/OM/CM and F01 retain their scopes.

