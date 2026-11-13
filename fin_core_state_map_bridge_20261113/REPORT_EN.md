# A concrete A11-to-RP state-map bridge

## Verdict

The A11 canonical 12+1 field template admits an exact two-field RP embedding
under explicit additional coefficient and kinetic assumptions. This supplies
a consistent state map, not a derivation of those assumptions from the strict
FIN kernel. A nonuniform Fourier eigenmode supplies a countercontrol: linear
spectral closure alone does not ensure nonlinear closure.

Status: CONDITIONAL_EXACT_EMBEDDING; STRICT_COEFFICIENT_SOURCE_OPEN.

This is the next bounded bridge calculation after the PDE campaign. It is
neither another audit of AW nor a repetition of the generic source inventory.
Its new result is a necessary-and-sufficient tangency test for a fixed linear
one-mode embedding, an exact positive example, and an explicit failed example.

## Actual candidate and provenance

The candidate is [A11](../fundamental_action_reconstruction/A11_NADSOLITON_SINGLE_KERNEL_CORE_LAGRANGIAN_AND_EMERGENCE_PACKET.md),
sections 4.1-4.3, checked against the potential constructed in
[QW-2165 source](../QW_2165_L13_EXHAUSTIVE_CANONICAL_EOM_GATE.py).
The [QW-2163 report](../RAPORT_QW2163_L13_FULL_CANONICAL_ACTION_VARIATIONAL_GATE.md)
explicitly leaves the complete-action all-orders theorem open. A11 likewise
states that the local coefficient families are not strict-derived from the
kernel alone. Its title and use of the word core do not remove that boundary.

K1/K2/F2/F3, S2, SUMMARY_GROK and live AGENTS were consulted. The current
state-map-first discipline takes precedence over older generic bridge routing.
No legacy physical-role transfer, selector closure or new lower informational
substrate is assumed. The nadsoliton remains the primordial informational
object. The pre-existing A11 Markdown wrapper and historical source files are
left unchanged. The old QW generator was read, not executed; this work does not
rerun its Lean wrapper or claim new formal verification of its contents.

The selected kernel has the form K_sg(d)=cos(omega*d+theta)/(1+beta*d^eta).
Here d is the internal ring index distance, not field amplitude or physical
time. Its nonlinear dependence on d does not, by itself, specify a cubic or
quintic field force. This study tests its actual place in the A11 action.

## State map and exact closure test

Write the supplied potential as

    V = mp Phi^2/2 + lambda Phi^4/4
        + Psi^T A Psi/2
        + sum_i [c4_i Psi_i^4/4 + c6_i Psi_i^6/6
                 + y_i Phi^2 Psi_i^2],
    A = diag(m_i) + (K + K^T)/2.

The symbols m_i and mp denote signed mass-squared coefficients. The directed
kernel contributes only its symmetric part to this real scalar potential;
there is no additional factor of two. This agrees with the bidirectional
coefficients in the QW-2165 Euler-Lagrange calculation.

Assume a common canonical wave kinetic term in one space dimension and a
constant real unit vector u. Define

    Phi = phi,   Psi_i = u_i psi,   P_perp = Id - u u^T.

The kinetic normalization is canonical because sum_i u_i^2=1. For arbitrary
phi and psi values, the full field force is tangent to this two-field
subspace if and only if all four vector conditions hold:

    P_perp A u = 0,
    P_perp diag(y_i) u = 0,
    P_perp (c4_i u_i^3)_i = 0,
    P_perp (c6_i u_i^5)_i = 0.

Proof: the transverse force is the polynomial

    psi P_perp A u + 2 phi^2 psi P_perp diag(y_i) u
    + psi^3 P_perp(c4_i u_i^3)_i
    + psi^5 P_perp(c6_i u_i^5)_i.

The four distinct monomials are independent. Their coefficient vectors must
vanish, and their vanishing is sufficient. Constant u also makes all
transverse kinetic derivatives zero. With tangent initial data, the embedded
solution solves the full equations for its interval of existence; uniqueness
within the chosen classical solution class gives invariance. This does not
prove global existence or transverse stability. The criterion concerns the
whole fixed linear subspace, not every special trajectory or nonlinear map.

When it closes, the exact restricted potential has coefficients

    a = u^T A u,
    g4_eff = sum_i c4_i u_i^4,
    g6_eff = sum_i c6_i u_i^6,
    y_eff = sum_i y_i u_i^2.

Merely substituting a mode into the action gives these projections even if
closure fails. In that case the result is a Galerkin model, not an exact
subsystem of the full equations.

## Positive example: uniform ring mode

Take u_i=1/sqrt(12), the circulant ring kernel of A11, and equal local
coefficients m_i=m, c4_i=c4, c6_i=c6, y_i=y. These equalities are added
premises; the general A11 template allows independent families.
Let kappa=sum_(d=1)^11 K_sg(d). Then all four closure conditions hold and

    V_eff = mp phi^2/2 + lambda phi^4/4
            + (m+kappa) psi^2/2
            + c4 psi^4/48 + c6 psi^6/864 + y phi^2 psi^2.

The RP potential used in the simulations, up to an irrelevant constant, is

    V_RP = (phi^2-1)^2/2 + (w-2+6 phi^2) psi^2/2
           + psi^4/4 + psi^6/6.

It is reproduced exactly by the following additional choices:

| A11 template datum | Required choice | Status |
|---|---|---|
| mp | -2 | Target matching, not sourced |
| lambda | 2 | Target matching, not sourced |
| y | 3 | Target matching, not sourced |
| c4 | 12 | Target matching, not sourced |
| c6 | 144 | Target matching, not sourced |
| m | w-2-kappa | Target matching, not sourced |
| Psi preparation | Uniform normalized mode | Supplied sector preparation |
| Kinetic operator | Common canonical wave operator | Supplied dynamical premise |

In particular V_(psi^4)=6, V_(psi^6)=120 and V_(phi^2 psi^2)=12 at psi=0,
where the subscripts denote differentiation orders. For general RP couplings,
the same map gives c4=12 g4_RP and c6=144 g6_RP.

The full 13-component gradient identity is checked symbolically, not only
the reduced action. Thus the existing PDE trajectories have a conditional
lift into this A11 specialization. Their existence in that specialization
does not show that strict FIN selects it. The uniform mode does not select
an axis in the degenerate pair1 sector and does not discharge QW-2191.

## Negative control: a linear eigenmode leaks nonlinearly

Let v_i=cos(2*pi*i/12)/sqrt(6). It is a unit eigenvector of every symmetric
circulant ring matrix. Uniform local coefficients therefore close the linear
part, but the cubic force has the exact nonzero transverse component

    ||(Id-vv^T) v^(componentwise 3)||^2 = 1/576.

For c4 != 0 the cubic transverse force has magnitude |c4 psi^3|/24.
The fifth-order term cannot cancel that coefficient identically for all psi.
The projected numbers g4_eff=c4/8 and g6_eff=5 c6/288 therefore do not,
by themselves, define an exact two-field reduction on this mode. This is one
explicit obstruction, not a classification of all Fourier phases or sectors.

A second control changes only c4_0 to c4+delta in the uniform-mode example.
The squared transverse cubic coefficient becomes 11 delta^2/20736.
Thus the uniform reduction also needs an actual equality/source premise for
the local coefficient family; the ring kernel alone does not supply it.

## Why this does not fix the nonlinear coefficients

Keep the same kernel, homogeneous mass and coupling data, and change the
homogeneous quartic coefficient from 12 to 24. The uniform embedding remains
exact, the static psi=0 kink and both linearized kink operators are unchanged,
but g4_eff changes from 1 to 2. Both members remain positive-quartic examples.
The kernel's quadratic mixing term has zero direct fourth and sixth field
derivatives. This is a concrete A11 realization of the earlier general
spectral underdetermination witness, now including a consistent state map.

This does not rule out nonlinear effects from eliminating additional modes,
or a different strict-core dynamical construction. Such a route would need
an elimination/state-map theorem and its error bounds. No such theorem is
silently inferred here. Finite ring index distance is not the physical PDE
coordinate; the continuum, clock, locality and kinetic assumptions stay open.

## Verification and next bounded step

check_bridge.py checks exact SymPy identities, including all 12 carrier
forces, the order-field force, potential matching, derivative normalization
and the two transverse-leakage controls. results.json records every outcome.
This is executable symbolic verification plus the analytic tangency proof;
it is not independent peer review or a machine-checked global PDE theorem.
The input hash ledger binds the actual files used at the parent commit.

The next useful step is a transverse-stability test of this explicit
conditional uniform embedding: compute all omitted ring-mode fluctuation
operators about the RP kink and classify their negative/zero modes over the
existing w interval. An invariant subspace can still be unstable to omitted
perturbations. A negative transverse mode would invalidate interpreting the
embedded kink as a stable full-core excitation, even though the two-field
simulation is internally consistent. A nonnegative result would remain
conditional and would not source c4/c6. Do not present further coefficient
matching as their derivation.

AW stays conditional pending independent review. No strict kernel completion,
legacy-role transfer, physical units, empirical FIN confirmation, photons,
SM/GR or ToE closure follows from this package.
