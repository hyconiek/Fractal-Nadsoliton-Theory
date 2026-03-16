# F657 First Exported `S_sel_int` Strict‑Core Source Object Selector Output Operator Packet

Status: `F657_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_OUTPUT_OPERATOR_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `N548`, the seed-v1 lane exports an admissible strict-core selector
reduction operator on `pair1`:

```text
B_sel_s_sel_int_source_object_v1 -> R_sel_s_sel_int_source_object_v1.
```

The next honest downstream move is to export one explicit selector output
operator:

```text
R_sel_s_sel_int_source_object_v1 -> O_sel_s_sel_int_source_object_v1
```

without pretending that an emergent observer, strict-core selector closure, or
any global `QW-2191` discharge already exist.

## Definition

Let:

```text
R_sel_s_sel_int_source_object_v1 : pair1 -> Q_sel_v1
```

be the exported selector reduction operator with `Q_sel_v1 := span{q_+,q_-}`.

Define the selector-output carrier:

```text
Q_out_v1 := span{o_+, o_-}
```

and export the output operator:

```text
O_sel_s_sel_int_source_object_v1(q_+) := o_+
O_sel_s_sel_int_source_object_v1(q_-) := o_-.
```

## Hard limits (no false pass)

`F657` does **not** claim:

- emergent observer construction,
- admissible `S_sel_int`,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

