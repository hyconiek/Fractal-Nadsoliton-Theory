# P655 Current Exported `S_sel_int` Strict‑Core Source Object Selector Bridge Operator Probe

Status: `P655_EXECUTED_CURRENT_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_BRIDGE_OPERATOR_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Test whether the current repo exports one admissible strict-core selector bridge
operator on `pair1` derived from:

```text
E_orient_s_sel_int_source_object_v1 -> B_sel_s_sel_int_source_object_v1
```

without claiming `R_sel`, `O_sel`, selector closure, `QW-2191` discharge, or
ToE closure.

## Checked properties

1. the orientation input is admissible (`N546`),
2. `B_sel` is derived only from the exported orientation datum,
3. strict-core only,
4. symmetric and traceless on `pair1`,
5. exports an internal signed selector decomposition (`P_sel_plus/P_sel_minus`),
6. has a positive source-alignment witness (seed lane only),
7. bridge-ready for a later `R_sel`,
8. kernel-split-safe.

## Hard limits

This probe does not claim:

- admissible `S_sel_int`,
- actual `R_sel`,
- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

