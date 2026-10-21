# SP: a source-bound spatial cascade estimate on a growing interval

Date: 2026-10-21. Conditional same-author analytical proof, not a physical
source, independent review or complete nonlinear modulation theorem.

## 1. Actual operators, domains and norms

Keep the unchanged RP contract and prepared real/odd/even sector.
H=-partial_x^2+4-6sech^2 x, L=H+w, w in[5/4,7/4],
f=(sqrt3/2)sech^2 x, K=tanh x, F=K f^2. The odd continuous projection
Pc,o removes the shape g; the even projection Pc,e removes f. No
translation zero mode is inverted. Shape, continuum and both momenta
remain in the actual model and in CM's canonical map.

Let S_A(t)=sin(t sqrt(A))/sqrt(A). Define the actual two-propagator kernel

    T(a,b)=S_L(a) Pc,e [K f S_H(b) Pc,o F], a,b>=0.        (1)

For j>=0 use ||h||_(j,-5)=sum_(r<=j)||<x>^-5 partial_x^r h||_2, and

    N16(h)=sum_(r=0)^16 ||<x>^16 partial_x^r h||_1.

All fixed profiles below are Schwartz. N16 is a conservative sufficient
input norm, not an optimized regularity assertion for arbitrary nonlinear
radiation. The result will prove it on the localized input in (1);
it is not assumed for all unknown PDE remainders.

## 2. Exact mixed spectral matrix element

Use the paid Jost states
psi_k=e^(ikx)(3K^2-1-k^2-3ikK)/D_-(k),
D_-(k)=(1-ik)(2-ik), D_+(k)=conjugate(D_-(k)), dk/(2pi).
For real k,l the ordinary absolutely convergent integral is

    M(k,l)=<psi_k,K f psi_l>
     =-i(sqrt3/2)(k-l)(k+l)^2(k^2+l^2+4) S2(k-l)
                              /[16D_+(k)D_-(l)].        (2)

Here S2(q)=pi q/sinh(pi q/2), with S2(0)=2. To prove (2) set q=k-l
and I_n=integral K^n sech^2(x)e^(-iqx)dx. Integration by parts gives

    I0=S2, I1=-iq S2/2,
    I_(n+1)=[n I_(n-1)-iq I_n]/(n+2).

Expanding the two Jost polynomials leaves
9I5+(-6-3k^2-3l^2+9kl)I3+(1+k^2)(1+l^2)I1
+3iq[3I4+(kl-1)I2]. The recurrence reduces this EXACTLY to
-i q(k+l)^2(k^2+l^2+4)S2/16. No Fourier sampling is involved.

As a source-level consistency check, K f=(sqrt3/24)[partial_x,H] on
C_c^infinity, since H's potential derivative is12K sech^2 x.
Equation(2), not an unsupported manipulation of scattering boundary
terms, proves its continuous kernel. The free derivative can have a
delta contribution; it must not be dropped in an abstract commutator
argument. The multiplication kernel (2) itself is an ordinary function.

The matrix element vanishes at k=l and k=-l (double at the latter).
At k=0 it reproduces the conjugate J-transform in ET, J=r0 K f,
r0=(3K^2-1)/2. This checks its orientation and normalization.

In particular (1) has the actual double continuum representation

    integral dk dl/(2pi)^2 psi_k(x) M(k,l) M_F(l)
       sin(a sqrt(k^2+4+w))/sqrt(k^2+4+w)
       sin(b sqrt(l^2+4))/sqrt(l^2+4).

For this fixed smooth source the integrals, with the finite derivatives
used below, are absolutely convergent. Parities are inherited from F
and K f; no physical sector is enlarged by using the full real k,l axes.

These spatial zeros do NOT remove every temporal resonance. At the
already allowed w=3/2, l=1 choose k>0 with k^2=1+sqrt30. Then

    sqrt(k^2+4+w)=sqrt(l^2+4)+sqrt(w),

but M(k,l) and M_F(l) are nonzero. All real denominator factors are
nonzero and k>l>0. Thus a uniform inverse of this temporal phase on the
whole mixed coupling is not licensed by the factors in (2). This is a
mathematical countercontrol, not chosen physical particle kinematics.

This nonzero statement also survives the actual odd-input parity fold,
not merely an individual travelling Jost matrix element. Write
M_F(l)=i U(l)/D_+(l), U(l)>0 for l>0. The sum over l and -l in the
double representation has, apart from nonzero common factors, numerator

    (k^2-l^2)[(k+l)S2(k-l)-(k-l)S2(k+l)].

For k>l>0 the bracket equals pi(k^2-l^2) times
1/sinh(pi(k-l)/2)-1/sinh(pi(k+l)/2)>0. Hence it does not cancel on
the displayed temporal resonant curve. The even output standing-wave
fold changes its phase convention, not this nonzero spectral coefficient.

## 3. Local spatial estimates with their input norms paid

First, for U(b)=S_H(b)Pc,o F, for any fixed j and q=0,1,

    |partial_x^j partial_b^q U(b,x)|
                    <=C_j <x>^3(1+b)^(-3/2).             (3)

Here is a direct proof rather than an appeal to an unread dispersive
theorem. Write the actual continuum integral and change
k=y sqrt(4+y^2), so sqrt(k^2+4)=2+y^2. Its amplitude is Schwartz in y.
Derivatives in x multiply by powers of k or differentiate bounded K;
derivatives in y up to3 produce at most <x>^3. The fixed M_F is
Schwartz and vanishes at k=0. Thus the transformed amplitude B_x(0)=0.

Gaussian-regularized Fourier inversion gives, for b>0,

    integral B_x(y)e^(ib y^2)dy
      =e^(i pi/4)sqrt(pi/b)/(2pi)
                    integral Bhat_x(xi)e^(-i xi^2/(4b))dxi.

The zero B_x(0) removes the constant term. Use
|exp(-i xi^2/(4b))-1|<=xi^2/(4b) and
integral |xi|^2 |Bhat_x|<=C(||B_x''||_2+||B_x'''||_2).
This gives (3) for b>=1. The negative time-frequency exponential has
the conjugate estimate. For b<=1 the original spectral integral is
bounded. For q=1 insert the smooth frequency factor before the argument;
no derivative of a previously unspecified error is used.

Since K f and every derivative decay exponentially, (3) proves, for
h_b=K f U(b) and its b derivative,

    N16(h_b)+N16(partial_b h_b)<=C(1+b)^(-3/2).             (4)

Second, for real even Schwartz h and m=sqrt(4+w), a>=1,

    S_L(a)Pc,e h(x)
     = r0(x)<r0,h> sin(ma+pi/4)/sqrt(2pi m a)+R_a h(x),
    ||R_a h||_(2,-5)<=C a^(-3/2) N16(h).                 (5)

The time derivative has leading term
r0<r0,h>sqrt(m/(2pi a))cos(ma+pi/4) and a remainder bounded by
C a^(-3/2)N16(h) in ||.||_(1,-5). This is proved separately spectrally.

To pay uniformity in h and space, its distorted transform and its first
three k derivatives obey
|partial_k^r M_h(k)|<=C(1+|k|)^(-12)N16(h), r<=3.
Integrate by parts in x in the Jost formula: each derivative costs at
most12 spatial derivatives and3 moments of h, dominated by N16.
All bounded tanh coefficients and their derivatives are harmless.
Set k=y sqrt(2m+y^2). The resulting y amplitude and its first3
derivatives are L2 uniformly with bound C<x>^3 N16(h), including up to2
x derivatives and the one time derivative stated above. Subtract its
value at0 in the Gaussian identity just used. That value is precisely
r0(x)<r0,h>, with the coefficient in (5). The remainder is
C<x>^3 a^-3/2 N16(h); multiplication by <x>^-5 is L2 integrable.
For a<=1 the unexpanded spectral integral gives C N16(h).
This proves actual local spatial estimates, not only matrix elements.

## 4. The sharper joint kernel bound

ET paid the exact cross-kernel cancellation

    <r0,h_b>=<J,S_H(b)Pc,o F>=O((1+b)^(-5/2)),            (6)

also for its b derivative. Its actual distorted product starts at k^4,
not k^2: M_J/(ik^3)->-sqrt3/16, M_F/(ik)->1/4.
This uses the continuous projection; the shape remains separately.

Insert (4) and (6) into (5), and use the short-time bounds. For a,b>=0,

    ||T(a,b)||_(2,-5)+||partial_a T(a,b)||_(1,-5)
      <=C[(1+a)^(-1/2)(1+b)^(-5/2)
                      +(1+a)^(-3/2)(1+b)^(-3/2)].        (7)

The same inequality holds after one b derivative. This is the specified
spatial two-propagator estimate for the ACTUAL localized input F. It
does not assert a norm bound on arbitrary radiative forcing in energy
space or improve global radiating-field energy decay.

The slower outer factor cannot be deleted. In fact
q(b)=<J,S_H(b)Pc,o F> has q'(0)=-d^2/sqrt3<0, since <J,F>=0 and
<J,g><g,F>=d^2/sqrt3. For a fixed sufficiently small b>0, q(b)!=0.
Then (5) gives a nonzero a^-1/2 leading profile. On an infinite sequence
of oscillation peaks its weighted local norm is bounded below by
c_b a^-1/2. A uniform a^-3/2 bound for that fixed b is false.
This is a countercontrol on the continuous subkernel, not a claim that
a shape/continuum cancellation can be discarded from the full solution.

## 5. The actual amplitude bound, reused rather than postulated slow

For the prepared invariant sector, CM's actual relative Hamiltonian has
E2 coercive on X=(H1_odd x H1_even)x(L2_odd x L2_even).
The nonlinear energies satisfy |E-E2|<=C(rho^3+rho^4+rho^6),
rho=||state||_X, with C depending only on finite M4,M6 and the fixed band.
The actual initial energy is exactly

    E(0)=w epsilon^2/2+(9g4/70)epsilon^4+(4g6/77)epsilon^6.

Choose a fixed small ball where E>=c rho^2. For sufficiently small
epsilon, conservation and the boundary energy barrier prevent exiting
that ball, and give rho(t)<=C epsilon. This is the existing small-energy
stability mechanism, not a new independent source/stability campaign.
It applies to this exact invariant sector, not arbitrary25-field data.

For completeness the continuation payment is elementary here: polynomial
nonlinearities are locally Lipschitz H1->L2 in1D, the gapped wave group is
bounded in the energy norm, and the local fixed-point time depends on
that norm. The bound therefore gives all finite-time continuation.
Smooth initial data retain strong regularity on each finite interval by
differentiated energy/Gronwall estimates; no useful high-norm long-time
constant is inferred. In particular the actual old z=<f,psi> obeys

    |z(t)|<=C epsilon for every t>=0.                     (8)

No slow-envelope derivative bound or physical clock has been assumed.

## 6. A genuine growing-window bound for source components

Define causal analysis components using that ACTUAL z:

    V_z(t)=-6 integral_0^t S_H(t-r)Pc,o F z(r)^2 dr,
    U_cas,z(t)=72 integral_(0<=r<=s<=t)
                   T(t-s,s-r)z(s)z(r)^2 dr ds,
    U_4,z(t)=-g4 integral_0^t S_L(t-s)Pc,e(f^3)z(s)^3 ds.  (9)

The sign and factor72 are (-12)*(-6), as in the actual equations.
These are exact terms in the causal decomposition of those equations,
not new applied forces or independent waves. The leftover fields and
all leftover nonlinear terms remain and are not called higher order
without proof.

From (3), (7), (8) and the corresponding velocity estimates,

    ||(V_z,V_z,t)||_(H2_-5 x H1_-5)<=C epsilon^2,
    ||(U_cas,z,U_cas,z,t)||_(H2_-5 x H1_-5)
                              <=C epsilon^3[1+sqrt(1+t)],
    ||(U_4,z,U_4,z,t)||_(H2_-5 x H1_-5)
                              <=C |g4|epsilon^3 sqrt(1+t). (10)

For the double integral, change to a=t-s,b=s-r. The positive bound is

    integral_(a,b>=0,a+b<=t)
       [(1+a)^(-1/2)(1+b)^(-5/2)
             +(1+a)^(-3/2)(1+b)^(-3/2)] da db
       <=(4/3)[sqrt(1+t)-1]+4.

For the direct quartic term use (5)'s coarse a^-1/2 local bound with the
fixed profile Pc,e f^3. Differentiating in t introduces only the outer
cosine propagator: S_A(0)=0 removes the boundary term. No derivative of z
is needed. Therefore for every fixed S and t<=S epsilon^-2, all three
components in (10) are O_S(epsilon^2) in their stated local field norms.
Unlike OM's prescribed envelopes, this statement uses the actual amplitude
and requires no C2-slow-envelope assumption.

## 7. What remains unpaid

Equation(10) is NOT a bound on the whole actual continuum field u or v.
The undamped shape part, the exact eta remainder, continuum-continuum
products and all other nonlinear terms have not been closed by (10).
Their forcing does not automatically have bounded N16 norms. Nor does
a local weighted bound imply small global radiative energy or a whole-
trajectory approximation. The old NM error cannot simply be inserted
at epsilon^-2.

The next nonlinear bootstrap must retain those actual remainders,
CM's initial slip, temporal resonant curves and the precise local/global
norm distinction. Return the first exact unpaid estimate rather than
assume the remainder small. No averaged damping, all-time scattering,
natural action/clock/geometry, noise, physical source or F01 closure.

