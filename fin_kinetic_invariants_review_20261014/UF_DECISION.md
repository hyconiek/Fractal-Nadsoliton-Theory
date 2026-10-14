# Architect decision after KI: symmetric-start critical fluctuations

Date: 2026-10-14.
Decision: **SYMMETRIC_START_CRITICAL_FLUCTUATIONS**.
Status: research selection and proof route; UF execution has not started.

## Why move away from the localization fold now

The supplied three-rule fold mechanism has a controlled path limit and the
same-readout gate. Canonical specificity failed in the useful sense: the
normalized exit number belongs to the class, not uniquely to A7. Another
precision campaign or another normalized copy of that outcome is not P0.

The present fold theorem starts in a deliberately selected near-fold state.
The next intrinsic question is whether the already known symmetric parent
can produce a controlled non-Gaussian collective trajectory through its own
microscopic fluctuations, without seeding a particular daughter/sign.

Use the existing uniform state and its already known unique top Fourier mode
k6. Do not claim the pitchfork, its static barrier or its recursive descendants
as new. Reports131/137 explicitly reject unrestricted recursive splitting and
exact finite-N Markov parity; those negative results remain binding.

## Proposed derived limit, not a theorem yet

Let sigma_j=(-1)^j, m=sum_j sigma_j p_j, alpha=g Lambda6/12. On the existing
parity line `p_j=(1+m sigma_j)/12`, primary HB gives exactly

    m_dot=tanh(alpha m)-m.

At `g_c=12/Lambda6`, its cubic drift is `-m^3/3`. The other ten full occupation
modes are stable. A potential new **actual-generator** target is

    g_N=g_c+lambda N^(-1/2),
    Z_N=N^(1/4)m_N(N^(1/2)s),
    dZ=[(Lambda6/12)lambda Z-Z^3/3]ds+sqrt(2)dB.

The primary starts can be iid uniform labels (or the separately tagged
balanced-uniform type), giving Z_N(0)->0 and finite-copy fast moment control.
No direction/sign, auxiliary noise, drive or physical time is inserted.
Initialization and mathematical gain choice remain explicit premises.

The normalizable quartic candidate differs from the local unbounded cubic
fold potential. A proof still needs actual finite-N drift/noise, removal of
fast modes on the **longer** N^(1/2) horizon, tightness and localization.
The old N^(-1/4) fast tube no longer gives a vanishing exit bound using just
the previous second-moment argument. Use a correctly paid tube/estimate.

One particularly important gate is cancellation of a dangerous `m*fast`
soft drift term: its crude absolute scaling is O(1), not o(1). For HB, the
within-parity averages of every other feature vanish, so the soft derivative
in each fast mean is zero along the entire parity line. Pay that selection
rule explicitly; MET/BAR need their actual flux/regularity treatment.

## Do not confuse local dynamical equilibrium with global Gibbs equilibrium

The uniform instability occurs well above the first coexistence gain. At the
same fixed A, take the already admitted ordered coexistence state p_eq:

    F_g_c(p_eq)=-(g_c-g_eq) p_eq^T A p_eq/2<0=F_g_c(u).

This variational witness uses the known root and PSD operator, not an
assumption that one localized branch remains the only global phase at all g.
Hence a quartic limiting stationary law near u cannot automatically be the
full finite-N Gibbs limit. The order of N and long-time limits is part of the
programme and an important kill test for autonomous/persistent-bit claims.

## Selected chain

[English plan](../FIN_UNIFORM_CRITICAL_NEXT_CHAIN_20261014.md):

    UF-001 reuse/dedup uniform critical structure and initialization
      -> UF-002 actual soft generator and symmetry/fast-coupling gate
      -> UF-003 quartic critical path theorem on the correct horizon
      -> UF-004 stationary-versus-dynamical limit reconciliation
      -> UF-005 handoff and one architect checkpoint.

If the microscopic quartic candidate or fast elimination fails, report the
exact surviving term; do not repair it with another noise, memory/controller
or graph. If local long-time laws disagree with global Gibbs, retain the
metastable/local classification rather than claiming a generated universe.

Success would be a stronger theory of collective statistical emergence from
an unbiased reference preparation. It would not derive the initial state of
the universe, natural gain/kernel, space, matter propagation or quantum/gravity.
