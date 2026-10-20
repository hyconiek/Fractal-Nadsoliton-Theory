# DEL-003: singular fast collapse by an explicit variational argument

Date: 2026-10-20. This proves the needed slow resolvents/semigroups directly.
No dense-limit Mosco theorem is silently applied on the wrong Hilbert space.

## 1. Exact changing-space identification

Let nu(dz,dY)=varrho(z,Y)dzdY=rho(z)gamma_H(Y)dzdY. Use half-open fundamental
cells C_N(x) of the actual scaled type lattice, centred consistently at
xi_x=(Z_N(x),Y_N(x)). Their volume is v_N=12N^-23/4, soft width O(N^-3/4)
and fast diameter O(N^-1/2). On each fixed compact, every cell centre is
in actual D_N eventually. The actual uniform compact Stirling/partition
calculation in MT gives pi_N(x)/v_N->varrho(xi_x) uniformly there.

For f on D_N define on its cell, and zero off the union,

    U_N f(xi)=f(x) sqrt[pi_N(x)/(v_N varrho(xi))].            (1)

This is an EXACT isometry L2(pi_N)->L2(nu), not evaluation of an arbitrary
L2 representative at lattice atoms. Locally its density multiplier tends
uniformly to one. Let

    H_slow=L2(q_A dz) direct-sum L2(q_B dz),
    J f(z,Y)=f_(P(z,Y))(z).

J is isometric. For phasewise bounded smooth compact-z sections define
J_N f(x)=f_(P_N(x))(Z_N(x)). Static tightness, the uniform compact density
and the null interface show U_N J_N f->J f strongly. Extend through chosen
section sequences when needed; arbitrary null-set L2 values are never used.

The full-space limiting energy will be infinite off this closed slow range.
It is NOT densely defined on all L2(nu). That is why only its slow-range
variational consequence is used below.

## 2. Finite dual edge calibration and fast independence

Write c_N(x,y)=pi_N(x) A_N(x,y). For any antisymmetric edge field v,

    E_N(f)=1/2 sum_xy c_N(x,y)(f_y-f_x)^2,
    B_N(f,v)=1/2 sum_xy c_N(f_y-f_x)v_xy,
    Q_N(v)=1/2 sum_xy c_N v_xy^2,
    B_N=-sum_x pi_N(x)f_x sum_y A_N(x,y)v_xy,
    E_N(f)>=2B_N(f,v)-Q_N(v).                             (2)

All identities are exact, by reversibility and completing the square.
Suppose U_N f_N converges weakly to f and sup E_N(f_N)<infinity.
Take a smooth compact VECTOR field a(z,Y) supported strictly inside one
phase, and use v_xy=Delta Y dot(a_x+a_y)/2. All its relevant edges are
bulk opposite-parity edges for large N. By DEL-002 and Taylor expansion,

    (sum_y A_N v_xy)/sqrt(N)
         ->[div_Y a-(H Y).a]/24,
    Q_N(v)/sqrt(N)->(1/24)integral ||a||^2 dnu.             (3)

The coefficient fields converge strongly through (1) on their compact
support. Cauchy--Schwarz in (2) gives B_N/sqrt(N)=O(N^-1/4)->0. Passing
to the limit in its divergence identity therefore yields

    integral f[div_Y a-(H Y).a]dnu=0

for EVERY such vector field. Thus grad_Y f=0 distributionally in each
open phase. For almost every z, the A fibre is a disk times R8 and the B
fibre is the exterior of a disk in R2 times R8; BOTH ARE CONNECTED. Local
mollification/Fubini now give f=J f_slow almost everywhere. This step is
not merely a claim that a weakly harmonic function must be constant.

## 3. Sharp slow liminf, without assuming a boundary corrector

For a compact smooth SCALAR field a supported inside either phase take
v_xy=Delta Z(a_x+a_y)/2. DEL-002 gives

    sum_y A_N v_xy ->b(z)a+partial_z a,
    Q_N(v)->integral a^2 dnu.                             (4)

The mixed bracket is o(1). The Taylor remainder is O(N^-1/4) on support:
rates O(N^3/2), Delta Z=O(N^-3/4), Delta Y=O(N^-1/2).
Again use the exact divergence identity and weak convergence. Equation (2)
implies for all such a

    liminf E_N(f_N)
       >=-2integral (Jf_slow)(b a+partial_z a)dnu
                                              -integral a^2dnu.

Since partial_z varrho=b varrho, this functional is the interior weak
z derivative. Optimizing first over scalar multiples of a gives its L2
bound; the Riesz representation supplies that derivative. On each phase
it equals f'_P(z). Taking the supremum over compact fields exhausting the
two open phases gives the FULL weighted norm, including arbitrary approach
to the null boundary:

For clarity, existence of f'_P is not assumed here: both phases contain a
fixed fast test region for every z in a given compact interval (a small
ball for A, an annulus outside the maximum k3 radius for B). Product tests
there first identify the ordinary local weak derivative of f_P. General
interior tests identify the same derivative throughout the phase. Density
of compact tests in L2 of each open phase then gives the full dual norm.

    liminf E_N(f_N)>=E_slow(f_slow)
             =sum_P integral q_P(z)|f'_P(z)|^2 dz.          (5)

If that derivative is not square-integrable, the right side is infinity.
Boundary energy was not discarded as negative: (2) is a global lower bound
and interior tests exhaust the limit measure. No global fast spectral gap
or full-H coupling theorem is used.

## 4. Recovery with true phase leakage charged

For smooth bounded compact-z f_P choose f_N=J_N f. On retained SAME-phase
edges, Taylor expansion and G_N->2 in L1 give

    E_N(f_N) ->sum_P integral q_P |f'_P|^2.

The error from squared-increment Taylor terms is O(N^-3/4). Cross-phase
edges are confined to E_N; their order-one phase jumps cost at most
C_f N^(3/2)pi_N(E_N)->0. Thus they do not impose a spurious matching
condition f_A=f_B on the interface. Norm convergence is given by (1).

The limit domain is the direct sum of weighted H1 spaces with energy (5).
Compact smooth sections are form-dense: truncate at large |z|, charge their
weighted L2/gradient tails and cutoff derivative, then mollify on compact
intervals where q_P is smooth positive. A diagonal choice gives recovery
for every limit-domain function. No values of rough L2 functions at atoms
are assumed in this extension.

## 5. Resolvents and the precise semigroup consequence

For r>0 minimize E_N(u)+r||u||^2-2<g_N,u> with U_N g_N->Jg strongly.
Minimizers have bounded norm and energy. Any weak subsequential limit is
slow by (3); (5) and recovery identify its unique minimizer as
u=(r-L_slow)^(-1)g. At the finite minimizer,
E_N(u_N)+r||u_N||^2=<g_N,u_N>. The limit identity plus the two lower bounds
force norm convergence as well, hence U_N u_N->J u STRONGLY. No lost
mass in rare states or expanding tails is silently assumed away.

The limit generator on phase P is the conservative weighted-form operator

    L_P f=q_P^(-1)(q_P f')'
         =f''+[b+p'_P/p_P]f'.                              (6)

For each fixed t>0, approximate exp(t A_N) spectrally by
(I-t A_N/k)^(-k). The scalar approximation error is uniform over the
nonnegative spectrum and tends to zero: on compact spectral intervals use
the logarithm expansion; beyond a fixed R both tails are bounded by their
values at R, then let R increase. Resolvent convergence therefore composes
for each fixed k and yields

    U_N exp(t A_N)g_N -> J exp(t L_slow)g.                  (7)

This is a SLOW-RANGE strong semigroup theorem. It claims neither a C0
limit on all full fast functions, compact spectral convergence, every
eigenvalue, a numerical gap nor uniform-time error rates.
