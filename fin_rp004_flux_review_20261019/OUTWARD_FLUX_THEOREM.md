# RP004-C: positive original nonlinear outward flux on the frozen window

Date: 2026-10-19. **CONDITIONAL_OUTWARD_FLUX_PAID_WITH_REVIEW_PROOFS**.
An analytical classical-model theorem, not formal verification, independent
peer review or physical validation. Source C was not executed; this is new
dependent work after the proofs of A and repaired B in this review.

## 1. Model, quantifiers and prepared trajectory

Use exactly the supplied canonical action in the invariant real e1 sector:

    phi_tt-phi_xx+2phi(phi^2-1)+6phi psi^2=0,
    psi_tt-psi_xx+(a+6phi^2)psi+g4 psi^3+g6 psi^5=0,
    phi(0)=tanh x, phi_t(0)=0,
    psi(0)=epsilon sqrt(3/4)sech^2 x, psi_t(0)=0.

For every finite M4,M6 and `a in[-3/4,-1/4]`, `|g4|<=M4`, `|g6|<=M6`,
there is epsilon_*>0 such that for every `0<epsilon<epsilon_*` the smooth
solution satisfies the estimate below. The threshold can be chosen uniformly
in these specified bounds. This does not quantify a useful numerical
epsilon_* or physical time. Signed coefficients are not changed to favorable
values. The full 25-real-field parent contains this exact prepared trajectory:
all spectator components and imaginary parts remain zero by uniqueness.
NS orbital stability applies in its own stated full-parent positive-gap and
finite-coefficient regime. No uniformity over unbounded spectator nonlinear
coefficients is asserted by borrowing that separate theorem.

Freeze the original endpoints, not an optimized alternative:

    T=12log(1/epsilon), R=2log(1/epsilon)=T/6,
    F_exact(R,T)=integral_0^T [J(t,R)-J(t,-R)]dt,
    J=-phi_t phi_x-psi_t psi_x.

For each epsilon, R is constant for this whole time integral. No independent
wave packet, external drive, artificial pin, Phi6 term or transient
subtraction is permitted or used.

## 2. Theorem

Set `c0=625/1548288`. For all sufficiently small epsilon in the above class,

    F_exact(R,T) >= (c0/8) epsilon^4 T >0.                   (20)

This is transported canonical energy in a conditional 1+1 Cauchy model,
from an initially localized bound-mode preparation. It is not merely
global energy growth of a forced reference, a positive resolvent density,
or near-field oscillation relabelled as emission.

## 3. Proof and full error inventory

The [moving-cone proof](MOVING_CONE_THEOREM.md) gives the actual retarded
reference estimate `F_u(T/6,T)>=(c0/4)T` for T>=T0 uniformly in a.
It retains DC, switch-on, odd shape mode, both continuum roots and source
tails. Its final band is strictly inside the harmonic group cone; spatial
energy there has positive order T coefficient, not an assumed stationary
outgoing solution or flux-only numerical diagnostic.

The [corrected strong remainder proof](BOOTSTRAP_COMPLETION.md) yields
for the very same data and frozen T,R,

    F_exact >= epsilon^4(c0/4)T
             -C[epsilon^6(1+T)^5+epsilon^8(1+T)^7
                                      +epsilon^10(1+T)].     (21)

Divide the error by epsilon^4 T. Its three terms are
`C epsilon^2(1+T)^5/T`, `C epsilon^4(1+T)^7/T`, and
`C epsilon^6(1+T)/T`. All tend to zero. Also T tends to infinity and
`epsilon^2(1+T)^3` tends to zero. Choose a single epsilon_* small enough
to pay T>=max(6,T0), the strong bootstrap and error ratio<c0/8. This
establishes (20), including kink interference, bound-mode tails, all
nonlinear cross terms and both endpoint signs. There is no incompatible
order of limits: the uniform reference bound is inserted at this same T.

For an independent sign convention check, the original potential is

    V=(phi^2-1)^2/2+(a+6phi^2)psi^2/2
                           +g4 psi^4/4+g6 psi^6/6.

Its energy density `e=(phi_t^2+phi_x^2+psi_t^2+psi_x^2)/2+V`
satisfies `e_t+J_x=0` exactly by the two PDEs. Hence

    F_exact(R,T)=E_[-R,R](0)-E_[-R,R](T).

This agrees with the outward orientation of (17)--(21). Unlike the forced
reference, the exact nonlinear parent has no source-work term in this
identity. The order-epsilon^2 reference forcing is extracted from its own
internal nonlinear coupling; the expansion is paid on the stated finite
window. Transport order epsilon^4 log(1/epsilon) is much smaller than
initial excitation energy order epsilon^2. There is no perpetual prescribed
drive in the parent or inconsistency with total-energy conservation.

## 4. What is closed and what is not

This closes the original RP-004 **finite logarithmic-window** signed-flux
gate in the frozen sparse family. It does not prove all-time damping,
asymptotic stability, scattering, a radiation power law, absorption by a
second wall, active noncommuting-kernel sensitivity, physical photons or
the emergence of space/time/energy from the strict FIN ontology. The active
species is still blind to K23. F01/F03/F06/F07/F10 are unchanged.

The original source retains its correct historical A-GATE-open/C-waiting
status. The new proof is visibly located here and is not attributed to that
submission. Return to architecture as required by the task/strategy; no
automatic new field campaign is activated by (20).
