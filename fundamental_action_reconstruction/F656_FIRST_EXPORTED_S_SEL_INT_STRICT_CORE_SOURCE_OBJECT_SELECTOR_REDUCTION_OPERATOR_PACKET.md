# F656 First Exported `S_sel_int` Strict‑Core Source Object Selector Reduction Operator Packet

Status: `F656_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_REDUCTION_OPERATOR_PACKET_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `N547`, the seed-v1 lane now exports an admissible strict-core selector
bridge operator on `pair1`:

```text
E_orient_s_sel_int_source_object_v1 -> B_sel_s_sel_int_source_object_v1.
```

The next honest downstream move is to export one explicit selector reduction
map:

```text
B_sel_s_sel_int_source_object_v1 -> R_sel_s_sel_int_source_object_v1
```

without pretending that `O_sel`, strict-core selector closure, or any global
`QW-2191` discharge already exist.

## Definition

Let `E_orient_s_sel_int_source_object_v1 := (e_parallel, e_transverse)` be the
ordered orthonormal frame on `pair1 = span{c1,s1}`.

Define the selector-sector carrier:

```text
Q_sel_v1 := span{q_+, q_-}
```

and export the reduction map:

```text
R_sel_s_sel_int_source_object_v1(x)
  := (<e_parallel, x>) q_+ + (<e_transverse, x>) q_-
```

for `x` in the declared `pair1` carrier.

## Hard limits (no false pass)

`F656` does **not** claim:

- admissible `S_sel_int`,
- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

