# P177 Current Exported Preobserver Selector Bridge Operator Probe

Status: `CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_PREOBSERVER_SELECTOR_BRIDGE_OPERATOR_FROM_E_ORIENT_PRELM_V1_AFTER_P177`
As of: `2026-03-08`

## Goal

Test whether the current repo now exports one actual admissible strict-core
preobserver selector bridge operator:

```text
E_orient_preLM_v1 -> B_sel_preLM_v1
```

without pretending that `R_sel`, `O_sel`, selector closure, or `QW-2191`
discharge already exist.

## Result

The probe passes on current repo state.

The repo now exports:

```text
B_sel_preLM_v1 := |e_parallel><e_parallel| - |e_transverse><e_transverse|
```

on the typed basis `[u_T, u_L, u_M]`, with:

```text
P_sel_plus_v1  := (I_TL + B_sel_preLM_v1)/2
P_sel_minus_v1 := (I_TL - B_sel_preLM_v1)/2
```

and a positive source-side selector response:

```text
<pi_TL(S_preLM_strict_core_source_object_v1),
 B_sel_preLM_v1 pi_TL(S_preLM_strict_core_source_object_v1)> > 0
```

## What this proves

The current repo exports one actual preobserver selector bridge operator that

1. is derived only from the admissible orientation datum,
2. remains strict-core only,
3. remains preobserver only,
4. is symmetric,
5. is traceless on the topological-light plane,
6. exports an internal signed selector decomposition,
7. is bridge-ready for a later `R_sel`.

## Hard limits

This probe does not claim:

- actual `R_sel`,
- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
