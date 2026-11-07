# TB - threshold bands and a logarithmic direct-cubic bound

Date 2026-11-07. Parent d07d935136c6aba40b7a4067966156965eaf6c2e.
Self-audited conditional analysis in the supplied RP model.

## 1. A single spectral budget for the entire low-frequency output

Use FH's exact converted free phase forcing F_3=(F_1,F_2), including all
terms of CS1. Write nu(k)=sqrt(m^2+k^2) and define its energy density

    e_F(s,k)^2=nu(k)^2 |Fourier(F_1)(s,k)|^2
                                  +|Fourier(F_2)(s,k)|^2.

Use the unitary Fourier convention. Fix 0<epsilon<kappa<=1 and
T=S epsilon^-2. Define

    B_epsilon(F)=integral_0^T [
        epsilon^-1 integral_(|k|<=epsilon) e_F(s,k)^2 dk
        +integral_(epsilon<|k|<kappa) e_F(s,k)^2/|k| dk
                                      ]^(1/2) ds.     (TB1)

The retarded low-output trace obeys

    ||q_3,lo||L2(0,T)<=C_S B_epsilon(F_3).             (TB2)

To prove this, first fix an impulse X at source time s. For the core,
Cauchy-Schwarz in k and boundedness of the Fourier observation kernel
bound the squared trace integrated on any interval of length <=T by
C T epsilon ||P_core X||energy^2. This is at most
C_S epsilon^-1 ||P_core X||energy^2.

For the intermediate region use the two half-lines k>0 and k<0 and the
two signs of temporal frequency. On each branch the frequency nu is
one-to-one. Full-line temporal Plancherel bounds the squared homogeneous
observation by C integral e_X(k)^2/|k| dk in that region. The constant
is uniform in the inner cutoff epsilon. Summation of its disjoint dyadic
bands happens inside this spectral integral. Within a branch their
full-line temporal frequency supports are disjoint. The two spatial
half-lines can overlap in time frequency and are combined by the
square-of-two-summands inequality. Restriction to a finite observation
interval comes after this full-line bound.

Combine core and intermediate parts by another square-of-two-summands
inequality, then use Minkowski for the source-time integral. This proves
TB2. It preserves the full retarded source and its time dependence.
There is no assertion of temporal orthogonality between driven responses
on the finite observation window.

Equivalently, with disjoint bands I_j={2^j epsilon<|k|<=2^(j+1)epsilon}
truncated at kappa, the squared bracket in TB1 is comparable, within a
factor two, to

    epsilon^-1 ||P_core F||energy^2
                 +sum_j (2^j epsilon)^-1 ||P_j F||energy^2. (TB3)

Thus a square sum, followed by the source-time integral, is sufficient.
A separate triangle inequality over all bands would discard this gain.

## 2. Precisely what the current actual bound supplies

CS supplies sup_s ||F_3(s)||energy<=C epsilon^3. Its direct insertion
into TB1 gives B_epsilon<=C_S epsilon^1/2. The resulting squared trace
is O(epsilon), the admitted generic critical bound. FH pays |k|>=kappa
by O(epsilon^2) for each fixed kappa. Consequently a sufficient missing
estimate for the actual whole cubic gate is

    B_epsilon(F_3)<=C_S epsilon.                     (TB4)

Together with CS4 and FH10 it implies the desired whole local-continuum
integral bound O(epsilon^2). TB4 is an explicitly stated open obligation.
The criterion controls absolute spectral size and may be stronger than
a criterion exploiting cancellations between different source times.

For orientation, suppose on |k|<kappa one had the additional pointwise
spectral bound e_F(s,k)<=C epsilon^3 |k|^alpha, alpha>=0. Direct integration
in TB1 yields

    alpha=0    squared trace <=C_S epsilon^2[1+log(kappa/epsilon)],
    alpha>0    squared trace <=C_(S,alpha,kappa) epsilon^2. (TB5)

In the first case the core contributes a bounded constant inside the
spectral bracket and the intermediate integral contributes the logarithm.
In the second case the core is O(epsilon^(2alpha)) and the intermediate
integral is bounded by kappa^(2alpha)/alpha. Constants may diverge as
alpha tends to zero. These are sufficient upper bounds. No sharpness
claim for the logarithm is made.

## 3. An actual improved bound for the direct cubic component

The direct physical continuum-only cubic force in CS1 is

    G_D=-6v^2 u-g4 u^3.

The global phase bound gives ||u||H1+||v||H1<=C epsilon. In one dimension
Sobolev embedding and Holder give

    ||v^2 u||L1 <=||v||L2^2 ||u||Linfinity,
    ||u^3||L1 <=||u||L2^2 ||u||Linfinity,
    ||G_D||L1+||G_D||L2<=C_(g4) epsilon^3.             (TB6)

The exact conversion of the direct signed half-wave pair gives physical
phase source (0,P_e G_D). Its free Darboux conversion is
(0,Wcal-adjoint G_D), because Wcal-adjoint P_e=Wcal-adjoint.
The conversion concerns this direct component only. CS1's other cubic
terms remain in the full F_3 used in TB1.

We need boundedness of Wcal-adjoint on L1, which follows from its explicit
formula. Set d(k)=sqrt((1+k^2)(4+k^2)) and let M_l be Fourier multipliers
with symbols

    M_0=-k^2/d(k), M_1=ik/d(k), M_2=1/d(k).

Commuting one derivative past K gives exactly

    Wcal-adjoint G=M_0 G+3M_1(KG)+M_2((2-3sech^2)G). (TB7)

The symbols M_0+1, M_1, M_2 and their first k derivatives are all in L2.
They are smooth at zero and decay at infinity as O(k^-2), O(k^-1),
O(k^-2), respectively. If a symbol a and a' lie in L2, its convolution
kernel has finite L1 norm by Plancherel and Cauchy-Schwarz with weight
(1+x^2)^-1. M_0 consists of minus the identity plus such a kernel.
The bounded multiplication factors in TB7 therefore prove

    ||Wcal-adjoint G||L1<=C||G||L1.                  (TB8)

Apply the L1-to-Linfinity Fourier bound to TB6--TB8. The direct phase
source satisfies e_D(s,k)<=C epsilon^3 for all k. TB5 with alpha=0,
and FH's fixed-cutoff bound for its high output, give the actual component
estimate

    integral_0^T |q_D(t)|^2dt
                     <=C_S epsilon^2[1+log(1/epsilon)]. (TB9)

The constants include the fixed finite signed coefficient g4. This is an
improvement over the generic O(epsilon) component bound. The logarithmic
factor and the separate cascade components still require further control.

## 4. A concrete sufficient cancellation, with its unpaid assumptions

For a real integrable G, TB7 or the adjoint identity with Wcal 1=r0 gives

    Fourier(Wcal-adjoint G)(0)
                          =(2pi)^-1/2 integral r0(x)G(x) dx. (TB10)

Even parity alone does not set this integral to zero. If, in addition to
TB6, the direct force obeyed

    integral r0 G_D=0,
    ||(1+|x|)G_D||L1<=C epsilon^3,

then its free Fourier transform would be bounded by C epsilon^3 |k|.
To justify this, the multiplier kernels in TB7 have integrable first
moments. Apply the preceding kernel argument also to a' using a'' in L2.
The identity part and the bounded spatial factors preserve the first
moment. The Fourier difference from its zero value is bounded by
|k| times that moment. TB5 with alpha=1 would then give O(epsilon^2).
These two source-specific assumptions are open for the actual force.
They provide a possible sufficient route alongside temporal cancellation.

## 5. Research outcome and next task

TB1 is a complete dyadic/core budget on the exact converted CS cubic
forcing. TB9 pays a sharper actual direct-cubic estimate using global
energy and the explicit L1 Darboux map. The cascade terms require their
own endpoint control, and removal of the remaining logarithmic loss
requires additional spectral or temporal information.

Next audit the complete normal-form cascade in the low-output Fourier
norm used in TB1. Derive its endpoint mapping bounds with the K factors,
commutators and physical phase projection retained. Then evaluate the
combined threshold coefficient of the direct and cascade terms before
assuming a cancellation. Actual critical closure, modal-action estimates
and independent proof review remain open.
