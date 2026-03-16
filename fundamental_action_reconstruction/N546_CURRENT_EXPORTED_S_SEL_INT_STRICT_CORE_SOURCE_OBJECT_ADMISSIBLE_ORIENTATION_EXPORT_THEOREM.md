# N546 Current Exported `S_sel_int` Strict‑Core Source Object Admissible Orientation Export Theorem

Status: `N546_DISCHARGED_CURRENT_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_ADMISSIBLE_ORIENTATION_EXPORT_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-16`

## Theorem statement

Within the current repo state, the exported strict-core source object:

```text
S_sel_int_strict_core_source_object_v1
```

exports one admissible orientation datum:

```text
E_orient_s_sel_int_source_object_v1
```

on the `pair1` frame `span{c1, s1}`.

## Scope

This theorem is scoped only to:

- current repo state,
- `S_sel_int_strict_core_source_object_v1`,
- `E_orient_s_sel_int_source_object_v1`,
- the seven-clause orientation export contract from `F32`.

## Proof basis

`P654` discharges all seven clauses (as checked against the explicit `F654`
export), therefore `E_orient_s_sel_int_source_object_v1` is an admissible
strict-core orientation export from `S_sel_int_strict_core_source_object_v1`
on current repo state.

## Hard limits

This theorem does not claim:

- admissible `S_sel_int`,
- actual `B_sel`,
- actual `R_sel`,
- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

