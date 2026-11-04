# CN: exact continuum-source normal form with the phase projection retained

Date: 2026-11-04. Base: 670954107a2502bf88f277d526601a1a4579de86.
Conditional component analysis in the unchanged RP model. No full critical
trace/action theorem, formal verification or independent review.

## 1. The actual inputs and the LF target

Use the whole physical continuous fields r_phi=v, r_psi=u in
eta=cg+v, psi=zf+u. Let P_phi,P_psi be the odd/even continuous projections.
They commute with A_phi=H and A_psi=L, but NOT in general with the free
multipliers Omega_i=sqrt(-partial_x^2+m_i^2), m_phi=2, m_psi=sqrt(4+w).
Write V=6sech^2(x), so A_i=Omega_i^2-V. The fields and physical velocities
start at zero in both continuous sectors and remain in their ranges.

The exact projected equations are

    r_i,tt+Omega_i^2 r_i=V r_i+G_i+L_i+P_i C_i,        (CN1)
    G_phi=-6K(v^2+u^2),  G_psi=-12Kvu,
    L_i=P_i B_i-(1-P_i)G_i,
    B_phi=-6K(c^2g^2+z^2f^2+2cg v+2zf u),
    B_psi=-12K(cz gf+cg u+zf v),
    C_phi=-6eta psi^2-2eta^3,
    C_psi=-6eta^2 psi-g4 psi^3-g6 psi^5.

Thus G_i+L_i=P_i(G_i+B_i). In particular the finite-rank subtraction
of G_i is essential: the extended products are not themselves continuous.
Every term of B has a localized bound profile. Global small energy and
boundedness of the projections give, for every fixed N,

    ||r_i||H1+||r_i,t||2 <=C epsilon,
    ||L_i||L2_N<=C_N epsilon^2,
    ||C_i||2<=C epsilon^3.                             (CN2)

The projected C_i also obey the latter global bound. Weighted bounds on
the entire C_i are not assumed.

LF's even source is P_psi(G_psi+N_X), where

    N_X=-6v^2u-g4u^3-g6u^5.                            (CN3)

Its exact zero-data response u_X is the object being studied. All source
coefficients still use the whole actual u,v, not u_X alone.

## 2. Physical phase projection in free half-wave coordinates

For a physical phase (h,p), define J_i(h,p)^sigma=h+sigma i Omega_i^-1 p.
Its inverse is h=(h^++h^-)/2 and p=Omega_i(h^+-h^-)/(2i).
The actual continuum projection in these coordinates is

    (mathcal P_i h)^sigma
       =P_i(h^++h^-)/2
          +sigma Omega_i^-1 P_i Omega_i(h^+-h^-)/2.    (CN4)

This is J_i diag(P_i,P_i) J_i^-1. It is idempotent, bounded on paired H1
and their polynomially weighted versions, and commutes with the full
linear half-wave generator corresponding to A_i. It does NOT mean
applying P_i separately to each half-wave. In particular the half-waves
of an actual continuous physical state need not each belong to ran(P_i).

Boundedness follows from the weighted free multiplier estimates already
used in QP/VD and the finite-rank smooth profiles of P_i. For example,
for P=1-|f><f| the explicit commutator is

    [Omega^-1,P]h=f<Omega^-1 f,h>-(Omega^-1 f)<f,h>.

Here Omega^-1 f has every polynomial moment by the multiplier estimates.
There is no discarded commutator. Most usefully, for the forcing vector
I_i(F)^sigma=sigma i Omega_i^-1 F, one has exactly

    mathcal P_i I_i(F)=I_i(P_i F).                     (CN5)

The full generator commutation is obtained through J_i from [A_i,P_i]=0,
not by incorrectly asserting [Omega_i,P_i]=0.

## 3. Reuse the paid free quadratic symbol on the actual continuum inputs

Set h_i^sigma=J_i(r_i,r_i,t)^sigma. Their exact free equations have forcing

    (partial_t+sigma i Omega_i)h_i^sigma=F_i^sigma,
    F_i=F_V,i+F_G,i+F_L,i+F_C,i,
    F_V,i^sigma=sigma i Omega_i^-1 V r_i,
    F_G,i^sigma=sigma i Omega_i^-1 G_i,
    F_L,i^sigma=sigma i Omega_i^-1 L_i,
    F_C,i^sigma=sigma i Omega_i^-1 P_i C_i.             (CN6)

These have H1 sizes O(epsilon), O(epsilon^2), O(epsilon^2), O(epsilon^3),
respectively. F_V and F_L have arbitrary fixed polynomial output moments
with those same sizes. The input h has global H1 size O(epsilon).

For the even output use QP's cross-channel coefficient a=-3 and its
symbol b_(sigma;tau,rho)=1/(Omega_psi(p+q)Phi), with
Phi=sigma Omega_psi(p+q)-tau Omega_phi(p)-rho Omega_psi(q).
QP's proof on the mass band gives bounded symbol and derivatives, the
H1 bilinear estimate, and the one-localized-input weighted estimate.
No bound frequencies replace the free masses in this symbol.

Define a quadratic correction evaluated on the actual continuous phase:

    Q^sigma(h)=sigma a K sum_(tau,rho) B_b(h_phi^tau,h_psi^rho). (CN7)

This is an analysis correction for a source response, not a claim of a
new canonical coordinate system. Its paired H1 norm is O(epsilon^2),
and Q(0)=0 since the physical continuous phase is initially zero.
This does not reset QP's NONZERO whole-state odd initial shift: QP used
the whole eta,psi input, whereas CN7 uses v,u only.

## 4. Exact transformed response and the complete forcing ledger

Let y be the zero-data half-wave response to the UNPROJECTED source
G_psi+N_X under the full linear operator A_psi. It is only an auxiliary
linear response; no global small-energy estimate for y is assumed.
By CN5, mathcal P_psi y is exactly J_psi(u_X,u_X,t).
Set Z=y-Q(h). Write the full linear half-wave operator as

    (mathcal D z)^sigma=(partial_t+sigma i Omega_psi)z^sigma
                         -sigma i Omega_psi^-1 V z_pos.

The same algebra as QP, with the ACTUAL forcings CN6, yields exactly

    mathcal D Z=R_2loc+R_3loc+R_ge3,  Z(0)=0,           (CN8)

where each sum below is over tau,rho and uses the cross-channel b:

    R_2loc^sigma = sigma i Omega_psi^-1 V Q_pos
       +sum[ sigma i a [Omega_psi^-1,K](h_phi^tau h_psi^rho)
              -i a [Omega_psi,K]B_b(h_phi^tau,h_psi^rho)
              -sigma a K{B_b(F_V,phi^tau,h_psi^rho)
                         +B_b(h_phi^tau,F_V,psi^rho)} ],

    R_3loc^sigma = -sigma a K sum[
                 B_b(F_L,phi^tau,h_psi^rho)
                   +B_b(h_phi^tau,F_L,psi^rho)],

    R_ge3^sigma = sigma i Omega_psi^-1 N_X
       -sigma a K sum[
          B_b(F_G,phi^tau+F_C,phi^tau,h_psi^rho)
            +B_b(h_phi^tau,F_G,psi^rho+F_C,psi^rho)].    (CN9)

For signs and coefficients, the primitive identity is

    (partial_t+sigma i Omega_psi)B_b(h_phi^tau,h_psi^rho)
       =i Omega_psi^-1(h_phi^tau h_psi^rho)
          +B_b(F_phi^tau,h_psi^rho)+B_b(h_phi^tau,F_psi^rho).

Multiplication by K produces its commutator; subtracting Q from y
produces the plus sign on Omega^-1 V Q_pos. This accounts for each term
in CN9. No potential, projection source or pure-bound forcing is dropped.
The L_i term explicitly carries the bound profiles and the finite-rank
quadratic correction from CN1.

Apply mathcal P_psi to CN8, using its commutation with the FULL generator.
Then convert with J^-1, and apply the Darboux map to both physical phase
components. This gives a free continuous KG response with the projected,
converted CN9 forcing and zero initial phase. No claim that projected
half-wave momentum is a time derivative is used for the transformed Z.
The physical u_X velocity is reconstructed from Z+Q exactly.

## 5. What the estimates pay, and what they do not

QP's explicit commutator and one-localized-input bounds, together with
CN2/CN6, give for every fixed N

    ||R_2loc||_(H1_N pair)<=C_N epsilon^2,
    ||R_3loc||_(H1_N pair)<=C_N epsilon^3,
    ||R_ge3||_(H1 pair)<=C epsilon^3.                  (CN10)

The phase projection, J^-1 and Darboux maps preserve the corresponding
weighted/global H1-position x L2-momentum classes. For J^-1 this uses
Omega:H1_N->L2_N; it does not require an extra source derivative.
The bounds concern sources, not the global size of each retarded response.

Define q_corr as the compact observation of mathcal P_psi Q after J^-1
and Darboux, and q_2loc,q_3loc,q_ge3 as the observations of the separate
zero-data continuous responses to CN9. Linearity gives exactly

    q_X=q_corr+q_2loc+q_3loc+q_ge3.                     (CN11)

The correction obeys |q_corr|<=Cepsilon^2. LE applied to R_3loc gives
local phase size Cepsilon^3 sqrt(1+t). Hence at T=S epsilon^-2,

    integral_0^T (|q_corr|^2+|q_3loc|^2)dt<=C_S epsilon^2. (CN12)

This pays the correction and the history generated by localized quadratic
forcing inside the normal form, including the retained bound sources.
CT's free observation estimate separately gives the available generic
bound integral |q_ge3|^2<=C_S epsilon, which is TOO LARGE for closure.
Absolute LE applied to R_2loc gives only an O_S(1) trace integral.

Consequently the remaining exact reduction is

    integral_0^T |q_X|^2dt
        <=C_S epsilon^2+4 integral_0^T(|q_2loc|^2+|q_ge3|^2)dt. (CN13)

The extended quadratic term has been removed from the transformed source,
but two unpaid histories remain. This is not a completed bound on q_X.
There is no general theorem that a cubic source is negligible on this
critical interval. BT's cubic internal threshold is not a zero of QP's
free quadratic symbol and remains relevant to later source-specific work.

## 6. Next narrowly stated tasks

First prove or refute a mixed estimate for the explicit R_2loc operators
of the form ||R_2loc||_(H1_N pair)<=Cepsilon E(t), with E the physical
local continuum norm in LF. QP's existing global-input localization bound
alone does not prove this stronger assertion. It requires control with a
negative-weight input through the bilinear and commutator operators.
Only after that payment may LF's absorption mechanism be reused here.

Second derive temporal/spatial cancellation or spacetime control for the
EXPLICIT R_ge3 in CN9, retaining both its direct cubic terms and quadratic
cascades. CT's norm-class countercontrol does not disprove such an estimate
for this source, but it excludes closure by its global size alone.

All three CD modal action errors remain separate. Preserve actual initial
data, signed coefficients, physical/QP/CM distinctions, strict/legacy and
source-selection boundaries. No full damping, scattering, empirical
validation or new FIN dynamics is claimed.
