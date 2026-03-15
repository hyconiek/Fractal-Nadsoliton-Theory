# P461 Current Strict Shannon Element-Order Reference `Z_n` Scope Extension Probe

Status: `P461_EXECUTED_ZN_SCOPE_EXTENSION_PROBE_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

The strict Shannon element-order reference lane currently exports a full `n=12` mode-index assignment by cutting each degenerate Fourier pair
via the defect:

```text
F_{2m}(ord) := Σ_{x=0}^{n-1} ord_{Z_n}(x) * exp(i * 4π m x / n).
```

For `n=12`, this defect is nonzero for all `m=1..5`, enabling the lane-scoped `O(2) -> Z2` cuts (`N480/N488/N496`, executed by `F454`).

This probe performs a **computational scope-extension check**:

- for a list of alternative carriers `n ≠ 12`,
- compute `ord_{Z_n}(x)` and the corresponding `F_{2m}(ord)` defects,
- report which pairs would be cut (`F_{2m} ≠ 0`) and which would remain degenerate (`F_{2m}=0`),

without promoting any generalization to theorem level.

## Inputs

- Pure group-structure computation on `Z_n` (no physics role-transfer claim).

## Output artifacts

- `fundamental_action_reconstruction/generated/p461_current_strict_shannon_element_order_reference_zn_scope_extension_probe.json`
- `fundamental_action_reconstruction/generated/p461_current_strict_shannon_element_order_reference_zn_scope_extension_probe_summary.json`

## Hard limits (no false pass)

This probe does **not** claim:

1. any theorem-level uniqueness / `O(2)->Z2` cut generalization beyond what is explicitly proven in the repo for the physical `n=12` lanes
   (note: the direction-freeness lemma `ord_{Z_n}` is `Aut(Z_n)`‑invariant is now recorded generally as `N503`, but this does not by itself
   promote any `n≠12` minimizer/assignment into theorem level),
2. any strict export of a mode-index assignment for `n ≠ 12`,
3. any discharge of `QW-2191` beyond the declared `n=12` lanes,
4. any ToE closure.
