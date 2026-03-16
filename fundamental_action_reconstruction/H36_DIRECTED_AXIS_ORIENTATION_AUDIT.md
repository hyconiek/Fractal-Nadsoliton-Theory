# H36: Directed Axis Orientation Audit

**Date:** 2026-03-16  
**Status:** `PASS_DIRECTED_AXIS_ORIENTATION_PRESENT__PREMISE_BASED_T164`

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

The strict core **now also supports** (premise-based via `T164` + `T171` discharge):

- an explicit sign-sensitive directed observable on `pair1` (`S_dir_pair1_strict_v1`, `F474`),
- an explicit global directed selector state object on `C_v1` descending to the projective state (`SelectorState_global_C_v1_directed_strict_v1`, `F474/N524`),

therefore the strict core does contain a directed orientation selection mechanism in the declared scope.

The strict core still does **not** support:

- a theorem that `u_psi0_pair1` is a physically selected directed axis rather than a chart-dependent representative,
- any `Aut(Z_12)`-invariant way to obtain the directed lift “for free” (ruled out by `N462`; the fixing datum is premise-based),
- strict-core selector closure (`S_sel_int`) or global discharge of `QW-2191`.

## Result

`H36_B2 := strict core exports a directed/sign-sensitive orientation selection mechanism on pair1 in the declared premise-based scope (T164+T171), while psi0 alone still supplies only an undirected axis representative`

## Hard limits

- No theorem-level PASS.
- No full-closure PASS.
- No claim that `psi0` selects a physically directed axis.
- No claim that `u_psi0_pair1` and `-u_psi0_pair1` are physically distinguished by strict core.
- No claim that `QW-2191` is discharged.
