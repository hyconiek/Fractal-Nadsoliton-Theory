# H37: Sign-Distinction State Audit

**Date:** 2026-03-16  
**Status:** `PASS_PARTIAL_SIGN_TRACKED_CONVENTION_LAYER_PRESENT_BUT_NO_STRICT_PHYSICAL_SIGN_DISTINCTION_OBSERVABLE`

## Goal

Check whether the current strict core contains any object that identifies `u` and `-u` inside `pair1 = (c_1,s_1)` as physically distinct states rather than coordinate-level representatives of the same undirected axis.

## Inputs

- `H31_PSI0_TO_PAIR1_REDUCTION_AUDIT`
- `H33_PAIR1_SELECTOR_TARGET_JUSTIFICATION_AUDIT`
- `H34_BASIS_COVARIANCE_TARGET_INDEPENDENCE_AUDIT`
- `H35_PAIR1_AXIS_SELECTION_AUDIT`
- `H36_DIRECTED_AXIS_ORIENTATION_AUDIT`
- `F467/N511`: lane-scoped oriented transport (α mod 2π) exists as a **tracked convention layer** (sign-tracked vector section), not a physical sign datum.
- `F470/N516`: global projective selector state object exported on `C_v1` (projector/span; residual sign is gauge at state level).
- `N517` (+ `F471`): even ord-reference weights (`ord_Z12`, `r_ord`) cannot supply a sign-distinction scalar of the form `Σ_x w(x) u_1(x)` on the current exported `pair1` sine axis.

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

Update (2026-03-16): the repo now also exports:

- a lane-scoped **sign-tracked convention layer** (vector-level oriented transport) on `{pair1..pair5}` (`F467/N511`),
- a strict global **projective** selector state object on `C_v1` (`F470/N516`),

but both remain compatible with residual sign as gauge unless an explicit sign-sensitive physical observable is exported.

Moreover, `N517` records a strict obstruction for one tempting route: the current strict `ord`-reference family (`ord_Z12`, `r_ord`)
is even under reflection and therefore cannot distinguish `u_1` from `-u_1` on the current exported `pair1` sine axis via an observable of the form `Σ_x w(x)u_1(x)`.

## Result

`H37_B2 := strict core exports sign-tracked convention-layer oriented vectors and a global projective selector state object, but still contains no strict sign-sensitive physical state object or observable distinguishing u from -u on pair1`

## Hard limits

- No theorem-level PASS.
- No full-closure PASS.
- No claim that sign reversal inside `pair1` is physically meaningful in strict core.
- No claim that `u` and `-u` are physically distinct selector states.
- No claim that `QW-2191` is discharged.
