# QP: quadratic far-field reduction with genuinely localized remainders

Date: 2026-10-21. Conditional same-author analysis, not a new FIN postulate.
This map is NOT asserted canonical. It is an exact local state-coordinate
change for the supplied PDE, retaining both fields and velocities.

## 1. The same source in free half-wave coordinates

Use the exact RP model and preparation, eta=phi-tanh(x), psi real,
eta odd and psi even, w=a0+2 in[5/4,7/4], finite signed g4,g6.
Set K=tanh(x), V=6sech^2(x), masses M_phi=2, M_psi=sqrt(4+w),
and FREE Fourier multipliers Omega_i=sqrt(-partial_x^2+M_i^2).
Their full linear field operators are Omega_phi^2-V=H and
Omega_psi^2-V=L. No bound mode or potential has been removed.

The exact second-order equations become

    eta_tt+Omega_phi^2 eta=V eta-6K(eta^2+psi^2)+N_phi,
    psi_tt+Omega_psi^2 psi=V psi-12K eta psi+N_psi,
    N_phi=-6eta psi^2-2eta^3,
    N_psi=-6eta^2 psi-g4 psi^3-g6 psi^5.                 (1)

K itself is not called an H1 state. For sigma=+/-1 define

    u_i^sigma=h_i+sigma i Omega_i^-1 h_i,t,
    h_phi=eta, h_psi=psi, h_i=(u_i^++u_i^-)/2,
    L_i^sigma=partial_t+sigma i Omega_i.

Thus L_i^sigma u_i^sigma=F_i^sigma. Its quadratic part is

    sigma i sum_(j,k,a in C_i) a Omega_i^-1 K
                                  sum_(tau,rho)u_j^tau u_k^rho,

where C_phi={(phi,phi,-3/2),(psi,psi,-3/2)},
C_psi={(phi,psi,-3)}. Its linear part is
F_i,lin^sigma=sigma i Omega_i^-1 V h_i, and its higher part is
sigma i Omega_i^-1 N_i. The branch sums have four terms; no attempt
clock, noise or microscopic update is introduced.

## 2. Uniform absence of FREE quadratic phase zeros

For p,q real, xi=p+q and signs sigma,tau,rho define

    E0=sqrt(xi^2+m0^2), E1=sqrt(p^2+m1^2), E2=sqrt(q^2+m2^2),
    Phi=sigma E0-tau E1-rho E2, b=1/(E0 Phi).             (2)

Use precisely the source mass triples
(2;2,2), (2;M_psi,M_psi), (M_psi;2,M_psi).
They obey |m1-m2|<m0<m1+m2. The Minkowski inequalities, here elementary,

    E1+E2>=sqrt(xi^2+(m1+m2)^2),
    |E1-E2|<=sqrt(xi^2+(m1-m2)^2),

follow from E1E2-pq>=m1m2 and E1E2+pq>=m1m2, respectively.
For example squared differences from m1^2 m2^2 equal
(pE2-qE1)^2 or (pE2+qE1)^2, and the left factors are positive.

For equal input signs the only small phase is E0-(E1+E2), and

    E0 |Phi|>=m0(m1+m2-m0).

For opposite input signs,

    E0 |Phi|>=[m0^2-(m1-m2)^2]/2.

For the three actual triples these bounds, and the easier large-phase
case, give uniformly for ALL momenta and all8 sign triples

    |E0 Phi|>=2, hence |b|<=1/2.                         (3)

M_psi>2 suffices for the cross-channel lower bounds. This is an
all-frequency proof, not extrapolation of a finite frequency grid.

Every mixed p,q derivative of b is also uniformly bounded. Here are
the needed details. Put S=E1+E2, D=E1-E2,
A=S^2-E0^2=m1^2+m2^2-m0^2+2(E1E2-pq)>0,
B=E0^2-D^2=m0^2-m1^2-m2^2+2(E1E2+pq)>0.
If equal input signs oppose sigma, b=sigma/[E0(E0+S)].
If equal input signs equal sigma, b=-sigma(S+E0)/(E0 A).
For opposite signs,

    b=(sigma+tau D/E0)/B.                               (4)

To bound derivatives without losing high-frequency cancellations, write
p=m1 sinh(theta), q=m2 sinh(varphi). Then
E1E2+/-pq=m1m2 cosh(theta+/-varphi).
All positive-order derivatives in p,q are bounded by a constant times
this positive quantity, since inverse rapidity derivatives are bounded
for masses bounded away from0. A,B have fixed positive lower bounds;
therefore their derivatives are bounded by C A,C B, and their inverses
have corresponding relative derivative bounds. The positive factors
E0,S+E0 have the same property. Formula(3) bounds the equal-sign ratios
themselves. Finally |D|<=E0, while every positive-order derivative of D
is bounded, so all derivatives of D/E0 are bounded. Formula(4) proves
the claim. All constants are uniform on the frozen compact mass band.

These are FREE masses2 and sqrt(4+w), NOT the bound frequency sqrt(w).
SP's actual bound/continuum temporal resonant curves are not contradicted
or removed: the localized potential below retains those interactions.

## 3. Bilinear maps and spatial localization without input moments

Let B_b(h,g) have free Fourier transform
integral b(p,xi-p) hhat(p) ghat(xi-p)dp/(2pi). For integer s>=1,

    ||B_b(h,g)||Hs
       <=C_s[||h||Hs ||g||H1+||h||H1 ||g||Hs].           (5)

Proof: use |b|<=1/2, the elementary frequency-weight product inequality,
Young's L2-L1 convolution inequality and ||hhat||1<=C||h||H1.
No unproved Coifman-Meyer application or derivative loss is hidden.

Define ||h||_(s,N)=sum_(r<=N)||x^r h||Hs. The bounded symbol derivatives
prove a ONE-localized-input version of (5):

    ||B_b(h,g)||_(s,N)
       <=C_s,N[||h||_(s,N)||g||H1
                         +||h||_(1,N)||g||Hs].           (6)

Indeed differentiate the output Fourier transform. Transfer derivatives
on ghat(xi-p) by integration by parts in p. One derivative becomes a
p derivative of b plus a derivative of hhat. Iteration puts every moment
on h or a bounded symbol derivative. This does not require moments of g.
Interchange the inputs for the other version.

The two remaining exact operators have real localization. For alpha=+1
or -1, Omega=sqrt(-partial_x^2+m^2), m in the source mass band,

    [Omega^alpha,K]:Hs -> Hs,N,
    Omega^-1 V:Hs -> Hs,N                              (7)

are bounded for every fixed integer s,N>=0. No input moment is assumed.
To prove it, use K'=sech^2, whose Fourier transform is S2(q), a Schwartz
function. In the commutator the distributional Fourier K becomes

    (1/i) [(Omega(xi)^alpha-Omega(p)^alpha)/(xi-p)] S2(xi-p).

The apparent pole is removable. The divided difference equals the
integral of the smooth derivative of (r^2+m^2)^(alpha/2) along the segment
from p to xi. All of its frequency derivatives are uniformly bounded.
Output Fourier derivatives, which measure x moments, act on this kernel
rather than on the input at p. Each resulting kernel is majorized by an
L1 Schwartz function of xi-p. With frequency weights and Young's
inequality this proves the first map in (7). For the second map the
kernel is Omega(xi)^-1 6S2(xi-p), and the same argument applies.
The possible delta convention for Fourier K is killed by the commutator.
This pays localization, not just unweighted boundedness.
The commutators are initially defined on Schwartz functions and extended
by these bounds; separate unbounded terms need not each act on Hs.

## 4. Actual invertible quadratic coordinate change

For each channel and signs use b from(2) and define

    Q_i^sigma(u)=sigma sum_(j,k,a in C_i) a K
                    sum_(tau,rho)B_b(u_j^tau,u_k^rho),
    v_i^sigma=u_i^sigma-Q_i^sigma(u).                    (8)

By(5), Q is quadratic and locally Lipschitz in H1. For a sufficiently
small H1 ball, u=v+Q(u) is a contraction on the corresponding ball.
This gives a genuine local inverse, not a formal inverse series. It is
near identity, with correction O(||u||H1^2). Tame higher regularity
follows from(5); for smooth actual fields all identities below hold
classically. SP's conserved-energy bound makes this H1 chart available
along the actual sufficiently small prepared trajectory for every finite
time. It does not give uniform high-Sobolev or weighted norms on long times.

The map respects reality u^- = conjugate(u^+) and the source parities:
b_(-sigma;-tau,-rho)=-b_(sigma;tau,rho), and the symbol is invariant
under (p,q)->(-p,-q). K is odd. Therefore eta remains an odd real
coordinate and psi an even real coordinate. No mode, continuum, velocity
or signed coefficient is discarded.

This map is NOT proved symplectic and must not be called canonical.
The reconstructed position/phase momentum need not obey the original
position-velocity identification. Physical fields and velocities are
recovered through the actual inverse and the half-wave definitions.
CM is a different canonical chart of the same source: composition is
possible through the physical state, but its modal coefficients and slip
must not silently be copied into v coordinates.

## 5. Exact transformed equation and remainder classes

For an individual channel,

    L_i^sigma B_b(u_j^tau,u_k^rho)
       =i Omega_i^-1(u_j^tau u_k^rho)
         +B_b(F_j^tau,u_k^rho)+B_b(u_j^tau,F_k^rho).

Apply L to K times this expression. The exact commutator gives, after
subtracting(8),

    L_i^sigma v_i^sigma-sigma i Omega_i^-1 V v_i,pos
                                =R_loc,i^sigma+R_ge3,i^sigma, (9)
    v_i,pos=(v_i^++v_i^-)/2.

Here the QUADRATIC part is explicitly

    R_loc,i^sigma =
       sigma i Omega_i^-1 V Q_i,pos
       +sum_(channels,branches) [
          sigma i a [Omega_i^-1,K](u_j^tau u_k^rho)
          -i a [Omega_i,K]B_b(u_j^tau,u_k^rho)
          -sigma a K ( B_b(F_j,lin^tau,u_k^rho)
                       +B_b(u_j^tau,F_k,lin^rho)) ],

and the remaining part is

    R_ge3,i^sigma =
       sigma i Omega_i^-1 N_i
       -sigma sum_(channels,branches) a K [
          B_b(F_j,ge2^tau,u_k^rho)+B_b(u_j^tau,F_k,ge2^rho) ].

All F, Q and N here are evaluated on the OLD u reconstructed from v.
F_ge2=F-F_lin is the actual quadratic-and-higher source, not a fitted
correction. Equations(6)--(7) prove, for rho_s=sum||u_i^sigma||Hs,

    ||R_loc||_(s,N)<=C_s,N rho_1 rho_s,
    ||R_ge3||Hs<=C_s rho_1^2 rho_s, rho_1 sufficiently small. (10)

The constants include the finite signed g4,g6 bounds. In particular ALL
quadratic remainders in this representation are spatially localized,
with their moments derived from unweighted Sobolev norms. The formerly
extended K eta^2/K psi^2/K eta psi products have not been dropped:
they were exchanged for the explicit coordinate correction and the
localized terms in(9). Full linear H/L dynamics and its bound modes remain.

For integer s>=16, N>=17 the first estimate pays N16(R_loc) by
Cauchy-Schwarz from weighted L2 derivatives. It does not assert rho_16
has an already paid uniform epsilon^-2 bound. Cubic terms need not have
positive moments and cannot be declared negligible: a size epsilon^3
term can accumulate an order-epsilon effect over time epsilon^-2.

In contrast, at s=1 the actual global small-energy bound from SP already
gives, for every fixed N and every finite model time,

    ||R_loc(t)||_(1,N)<=C_N epsilon^2,
    ||R_ge3(t)||H1<=C epsilon^3.                         (10a)

No positive moment of the unknown input fields was assumed to obtain
(10a). This is a whole transformed-forcing classification, not a bound
on its accumulated radiation or a new full-field asymptotic theorem.

## 6. The actual prepared data, including a NONZERO chart shift

At the physical Cauchy data, u_phi^+=u_phi^-=0 and
u_psi^+=u_psi^-=epsilon f. Let

    P0=sum_(tau,rho) B_(+;tau,rho)^(2;M_psi,M_psi)(f,f).

Then (8) gives EXACTLY, not only as a Taylor expansion,

    v_phi^+(0)=v_phi^-(0)=(3/2)epsilon^2 K P0,
    v_psi^+(0)=v_psi^-(0)=epsilon f.                      (11)

The new phase momenta are zero, but its odd position has a nonzero
quadratic shift. This is NOT incoming radiation or a new preparation.
The inverse map returns eta=eta_t=0, psi=epsilon f, psi_t=0.
One must calculate any new bound/continuum initial projections from(11),
not keep the old zero homogeneous data in this new chart.

The shift is genuinely nonzero. P0 is smooth and integrable by(6).
At output frequency0, p and -p have the same free energy
E=sqrt(p^2+M_psi^2). Summing the four symbols gives

    sum b=2/(4-4E^2)+1/2=(E^2-2)/[2(E^2-1)]>=13/34.

Here E^2>=21/4 and ||f||2=1. Parseval therefore gives

    integral P0 dx = P0hat(0)
       =integral (E^2-2)/[2(E^2-1)] |fhat(p)|^2 dp/(2pi)
       >=13/34>0.                                      (12)

These are FREE Fourier transforms of f, not its zero Darboux continuum
transform. Multiplication by K cannot annihilate a nonzero smooth
function except on the measure-zero set x=0. Thus (11) cannot be reset
to zero without changing the physical Cauchy state.

## 7. Paid dependency and exact remaining obligation

This pays a state map that removes extended quadratic products while
deriving localization of every quadratic remainder without assuming
radiation has positive spatial moments. It is an exact source-PDE
identity on its small-energy trajectory, not a new action or source law.

The coupled shape/radiation/bootstrap is still not closed. It must now
control the genuinely localized quadratic forcing, unweighted cubic
forcing, strong norms and their time integration, with data(11).
The localized potential retains the actual bound/continuum resonances;
nothing in the free phase estimate cancels SP's temporal countercontrol.
No full growing-window approximation, actual slow C2 envelope, averaged
damping, all-time scattering or physical F01/ToE closure is exported.

