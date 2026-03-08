# P178 Current Exported Preobserver Selector Reduction Operator Probe

Status: `CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_PREOBSERVER_SELECTOR_REDUCTION_OPERATOR_FROM_B_SEL_PRELM_V1_AFTER_P178`
As of: `2026-03-08`

## Goal

Test whether the current repo now exports one actual admissible strict-core
preobserver selector reduction operator:

```text
B_sel_preLM_v1 -> R_sel_preLM_v1
```

without pretending that `O_sel`, downstream completion, selector closure, or
`QW-2191` discharge already exist.

## Result

The probe passes on current repo state.

The repo now exports:

```text
R_sel_preLM_v1 : V_topo ⊕ L_int ⊕ M_int -> Q_sel_v1
```

with ordered selector-sector basis `[q_+, q_-]`.

For the current source object:

```text
R_sel_preLM_v1(S_preLM_strict_core_source_object_v1)
  = r_plus_v1 q_+ + r_minus_v1 q_-
```

where:

```text
r_plus_v1  > 0
r_minus_v1 = 0
```

## What this proves

The current repo exports one actual preobserver selector reduction map that:

1. is derived only from the admissible orientation datum and selector bridge,
2. remains strict-core only,
3. remains preobserver only,
4. exports an internal selector sector,
5. gives a positive plus-channel source response,
6. gives a vanishing minus-channel source response,
7. is bridge-ready for a later `O_sel`,
8. is kernel-split-safe.

## Hard limits

This probe does not claim:

- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
