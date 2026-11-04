# CS: isolate the exact cubic gate and audit its far-field resonance

Date: 2026-11-04. Base: b202477930e310a61e06dc6aac7b8f5271652b45.
Conditional analytical continuation. Actual source estimates and frozen
far-field reference calculations are explicitly separated below.
No formal verification, independent review or full damping theorem.

## 1. The exact cubic/remainder split along the actual trajectory

Retain CN's notation, phase projection, bound-forcing ledger and actual
physical continuous half-waves h_phi,h_psi. For the even output a=-3
and all sums are over input signs tau,rho. CN9's remaining forcing is

    R_ge3^sigma=R_3^sigma+R_hi^sigma,

    R_3^sigma=sigma i Omega_psi^-1(-6v^2u-g4u^3)
       -sigma a K sum[ B_b(F_G,phi^tau,h_psi^rho)
                         +B_b(h_phi^tau,F_G,psi^rho) ],

    R_hi^sigma=sigma i Omega_psi^-1(-g6u^5)
       -sigma a K sum[ B_b(F_C,phi^tau,h_psi^rho)
                         +B_b(h_phi^tau,F_C,psi^rho) ]. (CS1)

Here F_G is exactly the continuous-field quadratic forcing in CN6,
while F_C contains the projected WHOLE actual cubic/quintic physical
forcing. R_3 is a homogeneous cubic functional of the actual continuum
phase. R_hi starts at degree four in the whole small physical state.
Neither is an autonomous source; all inputs are the actual solution.
No replacement of whole u,v by their controlled response pieces is made.

The already paid global H1 input bound and bilinear estimate imply

    ||R_3||_(H1 pair)<=Cepsilon^3,
    ||R_hi||_(H1 pair)<=Cepsilon^4.                    (CS2)

Indeed F_C=O(epsilon^3) multiplied bilinearly by h=O(epsilon) is
O(epsilon^4); the direct quintic term is O(epsilon^5). These are global
bounds with no assumed positive spatial moments. The constants depend
on the fixed finite signed coefficients and original parameter band.
The conjugated phase projection, J^-1 and Darboux map preserve the
corresponding global phase bounds. All responses start at zero.

## 2. The higher-order history already meets the trace budget

Let q_3 and q_hi be the free-frame compact observations of the separate
retarded responses after those exact conversions. Then q_ge3=q_3+q_hi.
CT's observation estimate for a global phase source of size epsilon^p
implies

    integral_0^T |q_p|^2dt<=Cepsilon^(2p) T^2 sqrt(1+T).

At T=S epsilon^-2, p=4 gives

    integral_0^T |q_hi|^2dt<=C_S epsilon^3,             (CS3)

and the direct quintic alone gives O_S(epsilon^5). Thus the higher-order
history is smaller than the required O(epsilon^2) budget. It must be
retained in the exact identity, but no further cancellation is needed
for this particular estimate. This is an ACTUAL component bound.

Combining CS3 with ML15 reduces the remaining local-continuum gate to

    integral_0^T E(t)^2dt
                         <=C_S[epsilon^2+integral_0^T |q_3(t)|^2dt]. (CS4)

The full compact trace obeys the same right-hand form. The current
generic bound on q_3 is only O_S(epsilon), preserving ML's partial
whole-continuum result. The missing power is now located entirely in
the explicit cubic functional R_3, including its quadratic cascades.

## 3. Frozen far-field reference: no assumed cubic nonresonance

This section is a SYMBOL/REFERENCE calculation, not the actual prepared
trajectory or a proved limit theorem for it. Freeze K=kappa=+1 or -1,
V=0 and omit finite-rank bound-profile projections in the translation-
invariant far-field reference. Keep masses M_phi=2 and m=sqrt(4+w),
w in [5/4,7/4]. No such change is made to CS1--CS4 or to the PDE data.
The spatially constant trial below is not a global finite-energy state.

Set the first-order odd-field reference to zero and take

    u_1=(b exp(-imt)+conj(b)exp(imt))/2,
    v_1=0.

This tests the self-interaction coefficient of the even field, not the
entire multichannel cubic symbol. The free cubic temporal phase for
output sign + and input signs (+,+,-) at zero momenta is

    Phi_3=m-m-m+m=0.                                  (CS5)

Therefore a uniform nonzero cubic phase gap is unavailable even though
QP's FREE QUADRATIC phases have the previously proved gap.

The quadratic v response obeys v_2,tt+4v_2=-6kappa u_1^2. Its selected
stationary harmonic coefficients are

    v_2,0=-(3kappa/4)|b|^2,
    v_2,2m=(3kappa/[8(m^2-1)]) b^2.                   (CS6)

The coefficient at -2m is the conjugate. These coefficients use the
nonzero denominators 4 and 4-4m^2. Finite-time switch-on terms are not
claimed absent: they are excluded only from this selected stationary
harmonic calculation, and actual retarded histories remain in CS1.

Substitution in -12kappa v_2 u_1-g4 u_1^3 gives the coefficient of
b|b|^2 exp(-imt) in the physical reference source as

    A_self(w,g4)=9/2-9/[4(m^2-1)]-3g4/8.              (CS7)

The direct -6v^2u term has no cubic contribution when v_1=0; the quintic
term also has none. The two quadratic-cascade pieces in CS6 are both
required. Keeping only the direct g4 coefficient would miss them.

## 4. Check the same coefficient directly in CN's normal form

At zero momenta, in this reference h_psi^+=b exp(-imt),
h_psi^-=conj(b)exp(imt), h_phi^tau=0. The cubic cascade in CS1 still
survives because F_G,phi^tau=-3tau i kappa u_1^2. With output sign +,
QP's cross symbol becomes

    b_(+;tau,rho)=1/[m(m-2tau-rho m)].

The two exact branch sums are

    sum_tau tau b_(+;tau,+)=-1/m,
    sum_tau tau b_(+;tau,-)=1/[m(m^2-1)].              (CS8)

The resonant coefficients of u_1^2 h_psi^+ and u_1^2 h_psi^- are 1/2
and 1/4, respectively. Since -sigma a=+3 for sigma=+, the cascade
contribution is

    -9i[-1/(2m)+1/(4m(m^2-1))]
       =(i/m)[9/2-9/(4(m^2-1))].

Adding the direct cubic term yields exactly (i/m) A_self, agreeing with
the physical reference computation. K enters twice, so kappa^2=1 and
both far-field signs give the same coefficient. This is an algebraic
consistency check of the normal-form source, not a claim that one can
set the actual odd continuum to zero during nonlinear evolution.

## 5. Cancellation, phase, and the precise limitation

The selected coefficient vanishes only on the curve

    g4_star(w)=12-6/(w+3).                             (CS9)

At the band endpoints its values are 180/17 and 204/19. Its derivative
is 6/(w+3)^2>0, so a fixed g4 cannot cancel this reference coefficient
throughout the whole band. At w=4/3 the cancellation value is 138/13.
For g4=0 the coefficient is positive throughout the band. No cancellation
condition was imposed by the original signed-coefficient contract, and
no coefficient is retuned here. Vanishing of this ONE self-channel at
one parameter would not prove all mixed cubic channels vanish.

For the frozen resonant envelope equation alone,

    b'=(i/m)A_self |b|^2 b,

A_self is real, so d|b|^2/dt=0. The isolated resonant term is a phase
rotation in this reference equation, not a positive damping coefficient
or a proof of instability. It also does not show that the full actual
history has no amplitude change, radiation or accumulated trace.

A bounded translation-invariant cubic homological inverse obtained by
dividing this nonzero reference numerator by Phi_3 cannot cover the
self-resonant point generically. This rules out that particular blanket
nonresonance shortcut. It does not rule out resonant renormalization,
source-specific spacetime estimates or cancellation in the full actual
solution. The relation between this far-field reference symbol and an
actual distorted-wave threshold trace still has to be proved.

## 6. Next actual task

Treat the remaining R_3 as a structured resonant source rather than as
a uniformly invertible cubic phase. Derive its spectral/phase decomposition
with the actual potential, continuous projection and parity; identify
which resonant contributions affect phase and which can affect the
compact trace, and bound the remainders. A homogeneous reference phase
rotation alone is insufficient to import modified scattering or any
critical-time damping law into this model.

The three-channel modal action error remains separate. Preserve the
prepared physical trajectory, BT's internal mixed threshold at w=4/3,
whole-state QP initial shift, source/ontology/strict-legacy boundaries and
historical hashes. No physical-source or empirical FIN claim is added.
