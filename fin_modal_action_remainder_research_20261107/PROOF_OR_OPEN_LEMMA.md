# MA - canonical modal actions and a small integrated remainder

Date 2026-11-07. Parent 9cda9e97167259b104ed63df1c75a42ea1fb79fb.
Self-audited conditional analysis of the supplied RP model.

## 1. Fix the actual canonical chart and the paid estimates

Use CM's exact finite-rank map OLD=Phi3 composed Phi4(NEW), with new modal
coordinates (Z,P,C,R) and the unchanged physical continuum (v,v_m,u,u_m).
Put omega=sqrt(w), lambda=sqrt(3), and

    I=(P^2+w Z^2)/(2omega), J=(R^2+3C^2)/(2lambda).

These are observables on the actual PDE trajectory after CM reconstruction.
The transformed Hamiltonian is

    H_new=H2+H3m+Z4+Q4m+R5,
    Z4=A_z Jz^2+B_z Jz Jc+C_z Jc^2,
    Q4m=H4m+{H3m,G3}, ||dR5||<=C||state||^4.          (MA1)

CM proves these identities and uniformity on the parameter band. The
exact maps act as the identity on all continuum coordinates. The actual
global phase norm stays at most C epsilon. TP2 supplies, for T=S epsilon^-2,

    sup_(0<=t<=T) E(t)<=C epsilon,
    integral_0^T E(t)^2dt<=C_S epsilon^2,
    integral_0^T E(t)dt<=C_S.                         (MA2)

Here E is ML/LF's physical continuous local H1-position/L2-velocity norm
with weight <x>^-3, excluding the bound modes. The last inequality follows
by Cauchy-Schwarz and sqrt(T)=sqrt(S)/epsilon. The canonical map preserves
this norm exactly and the global smallness up to a uniform constant.

## 2. Retain the three radiative profiles and expand the other cubic terms

Use CD's profiles Fa=K f^2, Fb=K g^2 and G=K f g. The cubic interaction
containing exactly one continuum coordinate is

    H3rad=6 Z^2<Fa,v>+6 C^2<Fb,v>+12ZC<G,u>.          (MA3)

Their continuous projections can be inserted without altering the
pairings. This keeps the Z^2, C^2 and ZC channels with their exact factors.
The rest of H3m is

    H3two=6C<Kg,v^2+u^2>+12Z<Kf,vu>,
    H3pure=integral [2K v^3+6K v u^2] dx.             (MA4)

The last expression is independent of modal coordinates and has zero
Poisson bracket with I and J. With the canonical convention in CM,

    {I,H3rad}=-(12P/omega)[Z<Fa,v>+C<G,u>],
    {J,H3rad}=-(12R/lambda)[C<Fb,v>+Z<G,u>],
    {I,H3two}=-(12P/omega)<Kf,vu>,
    {J,H3two}=-(6R/lambda)<Kg,v^2+u^2>.               (MA5)

The momenta P,R here are canonical coordinates. No replacement by Z_t,C_t
is used. Those derivatives contain the full transformed interactions.

## 3. All nonleading action terms have integrated size O(epsilon^3)

Exponentially localized profiles pay all polynomial weights in MA2.
Thus the last two brackets in MA5 have absolute value at most
C epsilon E(t)^2, and

    integral_0^T (|{I,H3two}|+|{J,H3two}|)dt
                                                   <=C_S epsilon^3. (MA6)

The brackets with H2 and Z4 vanish exactly. Every modal-dependent monomial
of Q4m contains at least one continuum factor and at least one exponentially
localized profile. This follows for H4m by expanding the physical quartic
energy, and for {H3m,G3} because G3 differentiates only modal coordinates.
Pure-continuum terms have zero modal bracket. In each action bracket the
total degree remains four. With one continuum factor the bound is
C epsilon^3 E. With two it is C epsilon^2 E^2. With three, one remaining
modal factor gives C epsilon times a localized cubic continuum integral,
which is at most C epsilon^2 E^2 using global H1-to-Linfinity control.
Therefore

    |{I,Q4m}|+|{J,Q4m}|
                     <=C[epsilon^3 E+epsilon^2 E^2],
    integral_0^T (|{I,Q4m}|+|{J,Q4m}|)dt<=C_S epsilon^3. (MA7)

This treats the momentum-dependent quartic terms created by CM as well.
Finally |dI|+|dJ|<=C epsilon and the CM remainder gives

    |{I,R5}|+|{J,R5}|<=C epsilon^5,
    integral_0^T (|{I,R5}|+|{J,R5}|)dt<=C_S epsilon^3. (MA8)

The quintic force and all higher Taylor terms are included in R5. All
estimates are uniform for endpoints t<=T. Consequently the actual action
identities are

    I(t)-I(0)=W_I(t)+r_I(t),
    J(t)-J(0)=W_J(t)+r_J(t),
    W_I(t)=-(12/omega) integral_0^t P[Z<Fa,v>+C<G,u>]ds,
    W_J(t)=-(12/lambda) integral_0^t R[C<Fb,v>+Z<G,u>]ds,
    sup_(0<=t<=T)(|r_I(t)|+|r_J(t)|)<=C_S epsilon^3.  (MA9)

The remainder is o(epsilon^2) as epsilon tends to zero at fixed S. The
leading work has absolute integral at most C epsilon^2 integral E,
which is O_S(epsilon^2). Its sign and averaged constitutive law require
an estimate of the actual retarded continuum history.

## 4. Explicit physical observables and preparation

The corrected actions in physical variables are the exact pullbacks
I_corr=I composed (Phi3 composed Phi4)^-1 and likewise J_corr.
Their leading corrections are

    I_corr=I_phys-{I_phys,G3}+O(epsilon^4),
    J_corr=J_phys-{J_phys,G3}+O(epsilon^4).            (MA10)

Write G3=d[A z^2 r+B z p c+C_3 p^2 r+E_3 c^2 r+F_3 r^3],
with precisely CM's coefficients. The explicit cubic brackets are

    {I_phys,G3}=d[B omega z^2 c
             +(2omega C_3-2A/omega)z p r-(B/omega)p^2 c],
    {J_phys,G3}=d[lambda A z^2 c+lambda C_3 p^2 c
             +lambda E_3 c^3+(3lambda F_3-2E_3/lambda)c r^2
             -(B/lambda)z p r].                     (MA11)

The finite-rank map makes |I_corr-I_phys|+|J_corr-J_phys|<=C epsilon^3
uniformly on the actual small state. All initial shifts from CM are kept

    Z0=epsilon+O(epsilon^3), P0=R0=0,
    C0=-d A epsilon^2+O(epsilon^4),
    I(0)=omega epsilon^2/2+O(epsilon^4),
    J(0)=lambda d^2 A^2 epsilon^4/2+O(epsilon^6).      (MA12)

At w=3/2 the shape chart slip vanishes exactly as proved in CM. No
persistent C=O(epsilon^2) assumption is used for later times.

Replacing each canonical modal factor in W_I,W_J by its reconstructed
physical counterpart changes each product by O(epsilon^3). Pairing with
the continuum and integrating gives O(epsilon^3) by MA2. Thus MA9 also
has a physical-action version with the same remainder order and physical
momenta z_t,c_t. This statement follows from reconstruction, while the
canonical derivation always uses P,R as in MA5.

## 5. An integrated envelope estimate for the next transfer step

Define actual canonical interaction amplitudes

    a(t)=exp(i lambda t)(C+iR/lambda),
    b(t)=exp(i omega t)(Z+iP/omega).

The full Hamilton equations from MA1 give

    |a'|+|b'|<=C[epsilon E+E^2+epsilon^3].             (MA13)

The cubic mixed gradient has size epsilon E+E^2. The pure quartic
normal form contributes O(epsilon^3). Localized mixed quartic gradients
and the derivative remainder satisfy the displayed bound using global
smallness. This includes the difference between canonical momenta and
position derivatives. MA2 then yields

    ||a'||L2(0,T)+||b'||L2(0,T)<=C_S epsilon^2,
    ||(a^2)'||L2+||(b^2)'||L2+||(ab)'||L2
                                             <=C_S epsilon^3. (MA14)

Conjugate products have the same bound, so both mixed sum and difference
carriers are covered. These integrated bounds may be used in a retarded
response estimate. They do not impose a pointwise slow-envelope ansatz.

## 6. The exact remaining three-channel action gate

Let k_a,k_b,k_x be CD's reference coefficients, with their existing
normalization. To infer the candidate dissipative action law one still
needs, uniformly for t<=S epsilon^-2,

    W_I(t)+integral_0^t [k_a I(s)^2+k_x I(s)J(s)]ds=o(epsilon^2),
    W_J(t)+integral_0^t [k_b J(s)^2+k_x I(s)J(s)]ds=o(epsilon^2). (MA15)

MA9 pays the other Hamiltonian action remainders. MA15 concerns the
leading actual continuum work itself, whose bound is only O(epsilon^2)
from MA2. Its proof must construct the causal quadratic responses to
Z^2,C^2,ZC, retain source coefficients and switch-on transients, and
control their difference from the whole actual continuum in the work
pairings. CD's frozen densities alone do not supply MA15.

Next use MA14 to estimate that retarded-response transfer and identify
the first remaining mixed error. A critical damping law, independent
proof review and physical validation remain separate obligations.
