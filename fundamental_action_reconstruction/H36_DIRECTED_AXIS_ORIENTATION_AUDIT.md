# H36: Directed Axis Orientation Audit

**Date:** 2026-03-06  
**Status:** `PASS_PARTIAL_NO_STRICT_DIRECTED_AXIS_SELECTION`

## Goal

Check whether the current strict core contains any object that selects not only a coordinate-level axis inside `pair1 = (c_1,s_1)`, but also a directed orientation on that axis, i.e. distinguishes `u` from `-u`.

## Inputs

- `H30_KERNEL_INVARIANT_PSI0_ANCHOR_AUDIT`
- `H31_PSI0_TO_PAIR1_REDUCTION_AUDIT`
- `H33_PAIR1_SELECTOR_TARGET_JUSTIFICATION_AUDIT`
- `H34_BASIS_COVARIANCE_TARGET_INDEPENDENCE_AUDIT`
- `H35_PAIR1_AXIS_SELECTION_AUDIT`

## Audit

The strict core currently supports:

- a deterministic kernel-invariant anchor candidate `psi0`,
- a coordinate-level embedding `u_psi0_pair1 = cos(psi0)c_1 + sin(psi0)s_1`,
- a local mode chart `pair1 = (c_1,s_1)`.

The strict core still does **not** support:

- a theorem that `u_psi0_pair1` is a physically selected directed axis rather than a chart-dependent representative,
- a rule identifying `u_psi0_pair1` and `-u_psi0_pair1` as physically inequivalent,
- any strict-core sign-selection or directed-orientation object attached to `pair1`,
- any exported selector datum that chooses one orientation of the same axis over its opposite.

## Result

`H36_B1 := strict core supports only a coordinate-level undirected axis representative u_psi0_pair1 inside pair1 and contains no strict argument selecting a directed orientation on that axis`

## Hard limits

- No theorem-level PASS.
- No full-closure PASS.
- No claim that `psi0` selects a physically directed axis.
- No claim that `u_psi0_pair1` and `-u_psi0_pair1` are physically distinguished by strict core.
- No claim that `QW-2191` is discharged.
