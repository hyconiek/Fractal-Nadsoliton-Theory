# F655 First Exported `S_sel_int` Strict‑Core Source Object Selector Bridge Operator Packet

Status: `F655_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_BRIDGE_OPERATOR_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `N546`, the seed prefix:

```text
S_sel_int_strict_core_source_object_v1 -> E_orient_s_sel_int_source_object_v1
```

is now exported as an admissible strict-core orientation export (under `F32`).

The next honest downstream move is to export one explicit selector bridge
operator:

```text
E_orient_s_sel_int_source_object_v1 -> B_sel_s_sel_int_source_object_v1
```

without pretending that `R_sel`, `O_sel`, strict-core selector closure, or any
global `QW-2191` discharge already exist.

## Definition

Let:

```text
E_orient_s_sel_int_source_object_v1 := (e_parallel, e_transverse)
```

be the ordered orthonormal frame on `pair1 = span{c1,s1}` exported by `F654`.

Freeze the strict-core selector bridge operator as:

```text
B_sel_s_sel_int_source_object_v1 := |e_parallel><e_parallel| - |e_transverse><e_transverse|
```

acting on `pair1`, with the associated selector projectors:

```text
P_sel_plus_v1  := (I + B_sel)/2 = |e_parallel><e_parallel|
P_sel_minus_v1 := (I - B_sel)/2 = |e_transverse><e_transverse|.
```

## Hard limits (no false pass)

`F655` does **not** claim:

- admissible `S_sel_int`,
- actual `R_sel`,
- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

