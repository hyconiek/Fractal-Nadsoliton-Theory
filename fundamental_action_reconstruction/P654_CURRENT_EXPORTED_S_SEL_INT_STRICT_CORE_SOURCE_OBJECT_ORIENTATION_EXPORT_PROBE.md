# P654 Current Exported `S_sel_int` Strict‑Core Source Object Orientation Export Probe

Status: `P654_EXECUTED_CURRENT_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_ORIENTATION_EXPORT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Test whether the currently exported strict-core source object:

```text
S_sel_int_strict_core_source_object_v1
```

now supports one admissible orientation export:

```text
E_orient_s_sel_int_source_object_v1
```

under the seven-clause contract fixed in `F32`.

## Tested clauses

1. derived from the (exported) strict-core source object,
2. strict-core only,
3. internal orientation datum (or theorem-level equivalent),
4. selector-bearing without external anchor,
5. quotient / gauge safe,
6. bridge-ready for `B_sel`,
7. no silent kernel substitution.

## Hard limits

This probe does not claim:

- admissible `S_sel_int`,
- downstream `B_sel -> R_sel -> O_sel`,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

