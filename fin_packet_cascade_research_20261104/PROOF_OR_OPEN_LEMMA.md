# PC: a retarded odd-continuum bound for the full finite packet source

Date: 2026-11-04. Base: b75be24d1491e0766c23508a3e8d5531ab670567.
Conditional LINEAR packet-driven comparison result, not the prepared
nonlinear RP trajectory. Self-audited; no formal verification.

## 1. Target and result

Let P_R(t)=Re U_R(t) be TP's exact even continuous packet for L, of width
R>=1 and unit L2 position norm at t=0. For fixed S and 0<=t<=S R^2 define
Y_R as the zero-data odd continuous response

    Y_R,tt+H Y_R=-6 P_o[K P_R(t)^2],
    (Y_R,Y_R,t)(0)=0.                                (PC1)

P_o removes the shape mode. This component does not include or delete the
separate shape response of a full two-field calculation. Use the full
finite source in PC1, including its spatial tails; do not replace it by
TP's non-L2 leading profile or its negative-weight local source expansion.
We prove

    integral_0^(S R^2)||(Y_R,Y_R,t)||X_-6^2dt<=C_S,   (PC2)

with constants uniform in R and the original mass band. X_-6 is the
weighted H1-position/L2-velocity norm. The proof also isolates an exact
initial normal-form transient whose all-time squared local integral is
O(R^-1). That transient is retained, not reset to zero.

If the comparison input packet is multiplied by epsilon, the response
is multiplied by epsilon^2, so PC2 becomes C_S epsilon^4. Taking
R=epsilon^-1 gives a critical-time comparison estimate. It is not a bound
on the actual continuum generated from RP's zero-continuum preparation.

## 2. Packet bounds strong enough for the full spatial tail

Define the Fourier norm

    A_1(h)=integral (1+|k|)|hhat(k)|dk.

For the FREE half-waves of the physical packet,

    h^sigma=P_R+sigma i Omega_psi^-1 P_R,t,

TP's exact Fourier formula gives uniformly through t<=S R^2

    ||h^sigma||H1<=C_S,
    A_1(h^sigma)<=C_S R^-1/2,
    ||Omega_psi^-1 V P_R||H1_N<=C_(S,N) R^-1/2.       (PC3)

Here V=6sech^2 and N is any fixed nonnegative integer. These are global
Fourier norms and genuine positive-weight local forcing bounds; they
are not inferred from TP's negative-weight remainder.

Details of the Fourier bound: before applying the Darboux differential
operator, its inverse-multiplier packet is R^-1/2 G_(R,tau)(x/R), up to
the unit carrier exp(-imt), with tau=t/R^2. The rescaled Fourier amplitude
is a Gaussian times

    exp(-i tau q^2/[sqrt(m^2+q^2/R^2)+m])
                       /sqrt((1+q^2/R^2)(4+q^2/R^2)).

For tau in [0,S] its derivatives have uniformly bounded Schwartz
seminorms: each derivative of the displayed factors has a uniform
polynomial majorant, and the Gaussian supplies decay. The same holds
with the extra factor sqrt(m^2+q^2/R^2) for the physical time derivative.
Write the Darboux output as G''-3K G'+(3K^2-1)G, with derivatives in x.
The broad terms G and its derivatives have A_1=O(R^-1/2) or smaller by
scaling. Terms with a derivative of K, or K^2-1, are localized with
L2 and H1 sizes O(R^-1/2), also after the spatial derivatives needed here.

For K G' and its nonlocalized differentiated terms, use

    ||fhat||1<=C||f||2^(1/2)||f'||2^(1/2),

proved by splitting Fourier space at a positive frequency and optimizing
Cauchy-Schwarz. For K partial_x^j G, j>=1, the L2 norm is O(R^-j),
and its derivative has norm O(R^(-j-1)+R^(-j-1/2)). This yields a smaller
bound than R^-1/2. Apply the same argument to one more derivative for
A_1. All localized terms have the claimed bound directly. The free
Omega_psi^-1 multiplier preserves A_1, giving PC3 for both half-waves.
Global H1 follows from energy and the free half-wave map. The last bound
in PC3 follows from the packet's bounded local derivatives, exponential
V and the weighted free multiplier estimates, with no tail cutoff.

## 3. The exact quadratic correction and its localized forcing

Use QP's channel (phi;psi,psi), coefficient a=-3/2 and symbol

    b=1/[Omega_phi(p+q)
         (sigma Omega_phi(p+q)-tau Omega_psi(p)-rho Omega_psi(q))].

All needed symbol derivatives are uniformly bounded by the existing
free quadratic phase theorem. For this LINEAR input, the exact half-wave
forcing is only F_V^tau=tau i Omega_psi^-1 V P_R. There are no actual
nonlinear input forces to discard or approximate in this comparison.

Define, with the full sum over tau,rho,

    Q^sigma= sigma a K sum B_b(h^tau,h^rho).           (PC4)

Let y be the auxiliary unprojected full-H response to -6K P_R^2 in free
half-wave coordinates. As in CN, the exact identity is

    D_phi(y-Q)=R_loc,
    (y-Q)(0)=-Q(0),                                  (PC5)

where D_phi z=(partial_t+sigma i Omega_phi)z^sigma
                     -sigma i Omega_phi^-1 V z_pos and

    R_loc^sigma= sigma i Omega_phi^-1 V Q_pos
      +sum[ sigma i a [Omega_phi^-1,K](h^tau h^rho)
             -i a [Omega_phi,K]B_b(h^tau,h^rho)
             -sigma a K{B_b(F_V^tau,h^rho)
                              +B_b(h^tau,F_V^rho)} ]. (PC6)

All branch signs, potential terms and commutators are retained. Apply
CN's conjugated physical phase projection J_phi diag(P_o,P_o) J_phi^-1,
not componentwise P_o on the half-waves. It commutes with the full linear
generator. After converting to physical phase, PC5 becomes the exact
identity

    (Y_R,Y_R,t)(t)=Q_c(t)-U_H,c(t)Q_c(0)
                          +integral_0^t U_H,c(t-s)F_c(s)ds. (PC7)

Q_c is the projected physical phase corresponding to Q; its second
coordinate is not assumed to be the derivative of its first. F_c is the
converted projected PC6. The sum in PC7 reconstructs the physical velocity.
Since the packet initial field is nonzero, Q_c(0) is NOT set to zero.
This is a change in a comparison response decomposition, not QP's separate
whole-state chart or a modification of RP's preparation.

## 4. Bounds on the correction and the positive-weight source

Young's Fourier convolution estimate gives, for bounded symbol b,

    ||B_b(h,g)||H1<=C||h||H1 A_1(g),
    ||B_b(h,g)||W1,infinity<=C A_1(h)A_1(g).           (PC8)

In the second line W1,infinity means the function and its first derivative
are bounded. Each symbol derivative obeys the same estimates. The first
line uses (1+|p+q|)<=(1+|p|)(1+|q|), not a moment assumption.
Multiplication by K preserves H1 and bounded first derivatives. Thus

    ||Q_c(t)||X_0<=C_S R^-1/2,
    ||Q_c(t)||X_-6<=C_S R^-1.                         (PC9)

For the weighted momentum in PC9, use boundedness of Omega_phi from
H1_-6 to L2_-6 and the physical projection. No pointwise derivative of
the phase correction is required.

The commutator terms in PC6 also obey H1_N<=C R^-1. Indeed PC8 puts the
products and bilinear outputs in H1_-4 with size C R^-1, and ML8 maps
these through the explicit localized linear kernels into H1_N. This also
controls Omega_phi^-1 V Q_pos. For B_b(F_V,h), the one-localized-input
proof of QP can be sharpened to

    ||B_b(F_V,h)||H1_N<=C_N ||F_V||_(H1,N moments) A_1(h). (PC10)

Transfer all output-frequency derivatives to the localized input and
symbol, as in QP, and apply the first Young inequality in PC8. Only
bounded symbol derivatives and positive moments of F_V are used. PC3
therefore gives C R^-1. The symmetric input term is identical. Conversion
and projection preserve the result:

    ||F_c(t)||X_6<=C_S R^-1.                          (PC11)

This pays the positive weight needed by LE for the exact transformed
source. It does not misclassify the original broad source as localized.

## 5. Odd homogeneous local energy without weighted initial data

We also need to control Q_c(0) without assigning it an unpaid positive
weight. The supplied odd continuous operator obeys

    integral_R ||U_H,c(t)(a,b)||X_-6^2dt
                                          <=C||(a,b)||X_0^2. (PC12)

Here the integral is over the entire time line, and a,b are odd continuous
energy data. A direct proof uses the odd normalized generalized mode

    e_o(k,x)=[(3K^2-1-k^2)sin(kx)-3kK cos(kx)]
                              /sqrt((1+k^2)(4+k^2)), k>0.

For 0<k<=1, its weighted L2 and H1 norms are at most Ck, using
|sin(kx)|<=k|x|. For k>=1, those norms are at most C and C(1+k).
Expand the two time-frequency branches with nu(k)=sqrt(4+k^2), and
apply Plancherel in time after the monotone change k->nu. The Jacobian
is nu/k. Position contributes the weighted H1 mode norm, and velocity
contributes nu^2 times its weighted L2 norm. Their combination is bounded
by C[nu^2|ahat|^2+|bhat|^2] after multiplying by that Jacobian: the k^2
zero removes the threshold singularity; high frequencies cost no extra
energy derivative. Spatial integration uses the polynomial weight and
Fubini. Density gives PC12 for energy data. Constants of the half-line
spectral normalization are fixed and absorbed into C.

This argument applies to the ODD CONTINUOUS sector. It is not an even
threshold estimate, nor an assertion that the shape eigenmode decays.
It follows with PC9 that

    integral_0^infinity ||U_H,c(t)Q_c(0)||X_-6^2dt
                                                   <=C_S R^-1. (PC13)

## 6. Finish the finite-window bound and retain the outgoing history

LE's odd X_6->X_-6 decay is (1+t)^-3/2, integrable in time. PC11 gives

    sup_(0<=t<=S R^2)
      ||integral_0^t U_H,c(t-s)F_c(s)ds||X_-6<=C_S R^-1. (PC14)

Together with PC9 and PC13, square PC7 and integrate to get

    integral_0^(S R^2)||(Y_R,Y_R,t)||X_-6^2dt
                                     <=C_S[R^-1+R^-2(S R^2)]<=C_S.

This proves PC2. The nontransient two terms are locally O(R^-1);
the homogeneous transient is controlled in integrated local energy and
is not asserted to have that sharper pointwise bound. A coarse uniform
pointwise phase bound C_S R^-1/2 follows separately from global energy
for the transient and PC9/PC14.

No stationary outgoing wave is substituted for the response. PC7 retains
the whole causal integral, including DT's outgoing harmonic contribution,
finite switch-on and all spatial tails. The result is a bound, not a
leading-coefficient asymptotic comparison with DT's harmonic profiles.

## 7. Remaining bridge and numerical scope

Next estimate feedback of this retarded odd response into the even
packet/cubic history, retaining the extended spatial product and initial
transient. Then establish any link to the actual generated continuum
separately. The original nonlinear data contain no TP packet. PC2 alone
does not improve ML's O(epsilon) actual bound or close the O(epsilon^2)
cubic trace and three-channel action-error gates.

The accompanying numerical experiment evolves PC1 on a periodic Fourier
grid with centered time differences and the exact linear packet as source.
It checks time-step, grid and domain sensitivity and energy/source-work
balance. It is a finite numerical comparison, not the analytical proof,
an independent review or a simulation of the full nonlinear RP trajectory.
