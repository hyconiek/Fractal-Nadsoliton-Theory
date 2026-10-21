# ST: actual subcritical-time transfer, without a slow-envelope premise

Date: 2026-10-21. Conditional same-author analytical theorem.
Model, parameter band and Cauchy preparation remain exactly RP's.

## 1. Statements and norms

Set eta=phi-tanh x, omega=sqrt(w), w=a0+2 in[5/4,7/4],
f=sqrt(3/4)sech^2 x, g=sqrt(3/2)sech x tanh x.
The actual data are eta=eta_t=0, psi=epsilon f, psi_t=0.
Use the LOCAL phase norm X_-6 from LINEAR_WEIGHTED_PROPAGATION:
weighted H1 position and weighted L2 physical velocity. Constants below
are uniform in the source band and finite signed |g4|<=M4,|g6|<=M6.

For all sufficiently small epsilon and all finite t>=0, the actual
continuous components satisfy

    e_phi(t)=||(Pc,o eta,Pc,o eta_t)||X_-6
                                  <=C(epsilon^2+epsilon^3 t),
    e_psi(t)=||(Pc,e psi,Pc,e psi_t)||X_-6
                                  <=C(epsilon^2 sqrt(1+t)+epsilon^3 t). (1)

No hypothetical modulated forcing or independent wave is substituted.
In particular for every fixed0<alpha<2 and S<infinity,

    sup_(t<=S epsilon^-alpha) [e_phi(t)+e_psi(t)]/epsilon ->0. (2)

For the FULL prepared trajectory, including the actual shape and bound
mode coordinates, define

    Err(t)=||(eta,eta_t)||X_-6
       +||(psi-epsilon f cos(omega t),
                        psi_t+epsilon omega f sin(omega t))||X_-6.

Then

    Err(t)<=C[epsilon^2 sqrt(1+t)
                 +epsilon^3(1+t)^(3/2)+epsilon^4 t^2].   (3)

Consequently for every fixed0<alpha<4/3 and S<infinity,

    sup_(t<=S epsilon^-alpha) Err(t)/epsilon ->0.          (4)

The exponents are sufficient for THIS bound, not sharp dynamical
thresholds. Neither endpoint alpha=2 nor alpha=4/3 is admitted by these
little-o conclusions. No epsilon^-2 averaged damping law is claimed.

## 2. The actual global energy input and the two distinct charts

SP's small-energy barrier and continuation proof gives, on the exact
prepared invariant sector, ||physical state(t)||_(H1 x L2)<=C epsilon
for every finite t. It uses the true conserved relative Hamiltonian and
positive odd/even quadratic gaps; signed finite g4,g6 are retained.

QP is an exact generally noncanonical quadratic map in free half-waves.
Its proof pays, on that actual trajectory,

    ||R_loc(t)||_(H1,N)<=C_N epsilon^2 for every fixed N,
    ||R_ge3(t)||H1<=C epsilon^3,                         (5)

and its inverse correction in PHYSICAL phase energy is O(epsilon^2).
The latter follows directly from ||Q(u)||H1<=C epsilon^2 and the
bounded map Omega:H1->L2. No long-time strong Sobolev norm is used.

CM is a separate canonical finite-rank map, identity on all continuum
coordinates/momenta. Its energy jet and C1 remainder are used ONLY after
reconstructing the physical state. No CM coefficients are copied into
QP variables; no QP initial field is set to zero.

## 3. Low-regularity Duhamel in the QP chart

Write QP's variables v_i^sigma and define its phase coordinates

    A_i=(v_i^++v_i^-)/2,
    Pi_i=Omega_i(v_i^+-v_i^-)/(2i).

These Pi_i are not silently identified with A_i,t. QP's exact equation
instead gives

    A_i,t=Pi_i+F_i,pos,
    Pi_i,t=-(Omega_i^2-V)A_i+F_i,mom,
    F_pos=(R^++R^-)/2, F_mom=Omega_i(R^+-R^-)/(2i).       (6)

Split F into its actual localized and cubic parts. The operator Omega_i
maps weighted H1_N to weighted L2_N: commute each integer x moment using
bounded derivatives of its Fourier symbol. Equation(5) therefore gives

    ||F_loc||X_6<=C epsilon^2,
    ||F_ge3||X_0<=C epsilon^3.                          (7)

The localized source has phi odd and psi even. Apply the exact first-
order phase Duhamel formula for the FULL H/L wave group, not a time
derivative of a position Duhamel formula. This avoids any unproved
time derivative of R and uses precisely (7)'s regularity.

The new QP initial phi position is the actual O(epsilon^2) weighted
profile (3/2)epsilon^2 K P0, with zero phase momentum. Its continuous
part is therefore O(epsilon^2) in X_6. QP's initial psi is exactly
epsilon f with zero phase momentum; its even continuous part is zero.

The odd estimate from the linear theorem has integrable kernel:
integral_0^t(1+s)^(-3/2)ds<=2. The even estimate has integral at most
2sqrt(1+t). Use these on the localized part and global energy boundedness
on the cubic part. Recover the PHYSICAL state with the O(epsilon^2)
inverse correction. This proves (1). It neither assumes nor concludes
that the new QP phase momentum equals its position derivative.

For (2), the relative powers at t=S epsilon^-alpha are
epsilon, epsilon^(1-alpha/2), epsilon^(2-alpha), all vanishing for alpha<2.

## 4. Modal forcing estimate in the CM chart

Let a(t)=(Z,P,C,R) be CM's exact canonical bound coordinates. All have
size O(epsilon), and their continuum fields and momenta are precisely
the PHYSICAL continuous components just bounded in (1).
Let e(t)=e_phi(t)+e_psi(t).

The mixed cubic energy has the explicit form, with u,v denoting only
physical even/odd continuous POSITIONS in this paragraph,

    E3m=integral K[
       6C^2 g^2 v+6C g v^2+2v^3
       +6C g(2Z f u+u^2)+6v(Z^2 f^2+2Z f u+u^2)]dx.    (8)

Every bound-coordinate derivative brings a localized f or g. For example

    partial_Z E3m=12 integral K[C g f u+Z f^2 v+f v u],
    partial_C E3m=integral K[
       12C g^2 v+6g v^2+12Z g f u+6g u^2].

By1D Sobolev embedding, the full physical and CM positions have
unweighted H1 size O(epsilon). Since f,g and their products decay
exponentially, weighted Cauchy-Schwarz gives

    |grad_bound E3m|<=C epsilon e(t),
    |grad_bound^2 E3m|<=C e(t).                          (9)

This uses only local continuous position norms plus global O(epsilon)
H1, not strong derivative traces or radiation moments.
The same differentiation of the quartic mixed energy gives
|grad_bound E4m|<=C epsilon^2 e(t).
CM has Q4m=E4m+{E3m,G3}; its two differentiated factors, with
grad G3=O(epsilon^2), grad^2 G3=O(epsilon), and (9), therefore imply

    |grad_bound Q4m|<=C epsilon^2 e(t).                 (10)

Terms independent of bound coordinates have zero such derivative;
no unlocalized pure-continuum energy term has been estimated as local.

CM's exact energy is E2+E3m+Z4+Q4m+R5. Its pure-bound quartic Z4 has
gradient O(epsilon^3), and the actual energy-space C1 bound on R5 gives
bound gradient O(epsilon^4). The canonical bound equations consequently
are a free oscillator with uniformly bounded linear group plus forcing

    |forcing(t)|<=C[epsilon^3+epsilon e(t)+epsilon^4].   (11)

No mixed continuum interaction or canonical momentum was dropped.
The C1 remainder bound is uniform on the small ENERGY ball; it is not
the old fixed/short-window strong Taylor estimate from NM.

## 5. Integrating actual errors without a circular bootstrap

The actual CM initial coordinates differ from
(epsilon,0,0,0) by O(epsilon^2); their precise source slip is retained in
CM. Its exact inverse differs from identity by O(epsilon^2) in phase
energy, uniformly along the actual small-energy trajectory.

Apply finite-dimensional oscillator Duhamel to (11) and insert (1).
The modal error against the physical free bound preparation is at most

    C[epsilon^2+epsilon^3 t+epsilon integral_0^t e(s)ds]
       <=C[epsilon^2+epsilon^3(1+t)^(3/2)+epsilon^4 t^2]. (12)

Add the continuous-component errors (1) and use fixed profile f,g
weighted energy norms. This proves (3).
The argument is NOT circular: global conserved energy first gives (5);
the linear theorem then gives (1) for actual radiation; only then does
the independent CM modal estimate use those bounds. No assumed slowness,
small continuum bootstrap, or uniform rho_16 occurs.

Divide (3) by epsilon at t<=S epsilon^-alpha. The potentially dominant
relative powers are

    1-alpha/2, 2-3alpha/2, 3-2alpha,

all positive for alpha<4/3. This proves (4). The power4/3 comes from the
absolute time integration in (12), not a proved optimal physical timescale.

## 6. What the theorem changes, and what remains

This is an actual local transfer theorem for the unchanged nonlinear
prepared solution, not merely a given-envelope convolution lemma.
It covers the full local first-order field/velocity prediction on
subcritical times. All shape/continuum modes and both phase coordinates
were retained and their smallness derived. No new physical law was added.

The norm is weighted H1 x L2, NOT global H2 x H1. Thus this does not
backport a stronger trace/flux theorem into RP/NM. Nor does (4) reach
epsilon^-2, derive an averaged damping coefficient, all-time decay or
scattering. Those require the oscillatory/resonant time cancellations
lost in (12), and control of cubic accumulation. F01 remains unpaid:
the action, geometry, coefficients and preparation are still supplied.

