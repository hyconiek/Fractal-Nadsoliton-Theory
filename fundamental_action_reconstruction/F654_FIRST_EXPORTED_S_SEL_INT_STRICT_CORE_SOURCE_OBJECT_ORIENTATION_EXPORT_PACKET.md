# F654 First Exported `S_sel_int` Strict‑Core Source Object Orientation Export Packet

Status: `F654_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_ORIENTATION_EXPORT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `N545`, the minimal strict-core admissibility contract for:

```text
S_sel_int_strict_core_source_object_v1
```

is fully discharged (as a source-seed object only). The next honest strict move
is to export one explicit internal orientation datum:

```text
E_orient_from_S_sel_int_strict_core_source_object_v1
```

on the typed `pair1` target frame frozen by `F649`:

```text
E_orient_target_frame_v1 := (c1, s1).
```

## Output

`F654` exports one explicit orientation frame derived only from:

- `S_sel_int_strict_core_source_object_v1` (via `F647`), and
- its frozen typing target frame `(c1,s1)` (via `F649`),

as a strict-core internal export (no imported `psi0`, no external selector
control).

## Hard limits (no false pass)

`F654` does **not** claim:

- admissible `S_sel_int`,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

