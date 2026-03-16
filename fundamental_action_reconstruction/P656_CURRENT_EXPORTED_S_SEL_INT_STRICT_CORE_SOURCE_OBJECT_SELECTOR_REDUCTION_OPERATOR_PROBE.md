# P656 Current Exported `S_sel_int` Strict‑Core Source Object Selector Reduction Operator Probe

Status: `P656_EXECUTED_CURRENT_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_REDUCTION_OPERATOR_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Test whether the current repo exports one admissible strict-core selector
reduction operator on `pair1` derived from:

```text
B_sel_s_sel_int_source_object_v1 -> R_sel_s_sel_int_source_object_v1
```

without claiming `O_sel`, strict-core selector closure, `QW-2191` discharge, or
ToE closure.

## Checked properties

1. orientation input admissible (`N546`),
2. selector bridge operator admissible (`N547`),
3. reduction derived only from orientation + bridge,
4. strict-core only,
5. internal selector sector exported (`Q_sel_v1`),
6. positive plus-channel and vanishing minus-channel on the seed source response,
7. bridge-ready for a later `O_sel`,
8. kernel-split-safe.

## Hard limits

This probe does not claim:

- admissible `S_sel_int`,
- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

