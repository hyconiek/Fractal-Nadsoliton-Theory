# H37: Sign-Distinction State Audit

**Date:** 2026-03-06  
**Status:** `PASS_PARTIAL_NO_STRICT_SIGN_DISTINCTION_STATE_OBJECT`

## Goal

Check whether the current strict core contains any object that identifies `u` and `-u` inside `pair1 = (c_1,s_1)` as physically distinct states rather than coordinate-level representatives of the same undirected axis.

## Inputs

- `H31_PSI0_TO_PAIR1_REDUCTION_AUDIT`
- `H33_PAIR1_SELECTOR_TARGET_JUSTIFICATION_AUDIT`
- `H34_BASIS_COVARIANCE_TARGET_INDEPENDENCE_AUDIT`
- `H35_PAIR1_AXIS_SELECTION_AUDIT`
- `H36_DIRECTED_AXIS_ORIENTATION_AUDIT`

## Audit

The strict core currently supports:

- a legal coordinate-level representative `u_psi0_pair1`,
- a local chart `pair1 = (c_1,s_1)`,
- a deterministic angle candidate `psi0`.

The strict core still does **not** support:

- a sign-sensitive observable on `pair1`,
- a state object distinguishing `u` from `-u`,
- any parity-breaking selector datum attached to the `pair1` chart,
- any theorem that sign reversal on the local axis changes physical content rather than only coordinates.

## Result

`H37_B1 := strict core contains no sign-sensitive state object or observable on pair1 and therefore does not distinguish u from -u as physically different selector states`

## Hard limits

- No theorem-level PASS.
- No full-closure PASS.
- No claim that sign reversal inside `pair1` is physically meaningful in strict core.
- No claim that `u` and `-u` are physically distinct selector states.
- No claim that `QW-2191` is discharged.
