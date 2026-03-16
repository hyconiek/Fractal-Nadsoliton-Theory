# H38: Projective Selector State Audit

**Date:** 2026-03-16  
**Status:** `PASS_DIRECTED_CONTINUATION_SELECTED__GLOBAL_DIRECTED_STATE_OBJECT_EXPORTED__PROJECTIVE_STATE_RETAINED_AS_QUOTIENT`

## Goal

Check whether the current strict core supports the selector state on `pair1 = (c_1,s_1)` only at the projective/ray level, rather than as a physically individuated directed vector.

## Inputs

- `H31_PSI0_TO_PAIR1_REDUCTION_AUDIT`
- `H33_PAIR1_SELECTOR_TARGET_JUSTIFICATION_AUDIT`
- `H34_BASIS_COVARIANCE_TARGET_INDEPENDENCE_AUDIT`
- `H35_PAIR1_AXIS_SELECTION_AUDIT`
- `H36_DIRECTED_AXIS_ORIENTATION_AUDIT`
- `H37_SIGN_DISTINCTION_STATE_AUDIT`
- `F470/N516`: global projective selector state object exported on `C_v1` (projector/span semantics; residual sign is gauge at state level).
- `N501` / `N502`: residual `Z2` sign flips are gauge-irrelevant for the currently exported downstream span/projector objects (declared scope).
- `F474/N524`: exports a sign-sensitive directed observable on `pair1` and a global directed selector state object on `C_v1` (vector-level lift; `T171` discharged, premise-based via `T164`).
- `P475`: earlier decision packet selecting **projective-only continuation** (historical branch; superseded if directed continuation is selected).
- `P632`: professorial decision packet selecting **directed continuation** (treat the directed/vector-level state as physical in the declared scope).

## Audit

The strict core currently supports:

- a legal coordinate-level representative `u_psi0_pair1`,
- a local chart `pair1 = (c_1,s_1)`,
- a deterministic angle candidate `psi0`,
- a strict global **projective/ray-level** selector state object on `C_v1` (`F470/N516`),
- and (premise-based via `T164`) a strict global **directed/vector-level** selector state object on `C_v1` (`F474/N524`).

Taken together, these facts now support two explicit layers (no false pass):

1. **Projective layer (quotient):** the ray/projector state where `u` and `-u` are identified (still exported and still useful for sign-gauge-safe downstream objects).
2. **Directed layer (premise-based):** a vector-level state where `u` and `-u` are distinguished in the declared scope by an explicit sign-sensitive observable and an explicit fixing datum (`T164`).

The strict core still does **not** support:

- a strict-core selector closure claim (`S_sel_int`),
- a global discharge of `QW-2191`.

By `P632`, the current strict continuation proceeds under an explicit **directed** interpretation in the declared scope:

- treat the selector state as a directed/vector-level object (`F474/N524`),
- keep the projective state as the quotient shadow where appropriate,
- keep all generator/sign dependence explicit (premise-based via `T164`; no Aut-invariant canonicity claim by `N462`).

## Result

`H38_B2 := strict core exports both a global projective selector state object (quotient/ray layer) and a global directed selector state object (vector layer, premise-based via T164), with the directed layer descending to the projective one`

## Hard limits

- No theorem-level PASS.
- No full-closure PASS.
- No claim that the projective/ray-level state is already the physically correct selector state.
- No claim that `QW-2191` is discharged.
