# F91 First Exported Preobserver Selector Output Operator Packet

Status: `F91_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SELECTOR_OUTPUT_OPERATOR_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N198`, the next honest constructive move is:

```text
R_sel_preLM_v1 -> O_sel_preLM_v1
```

That is, export one actual preobserver selector-output map from the already
constructed selector reduction operator, without pretending that an emergent
observer, selector closure, or `QW-2191` discharge already exist.

## Reused selector reduction input

Reuse:

```text
R_sel_preLM_v1 : V_topo ⊕ L_int ⊕ M_int -> Q_sel_v1
```

with selector-sector basis:

```text
Q_sel_v1 := span{q_+, q_-}
```

and current source-side selector response:

```text
R_sel_preLM_v1(S_preLM_strict_core_source_object_v1)
  = r_plus_v1 q_+ + r_minus_v1 q_-
```

where:

```text
r_plus_v1  > 0
r_minus_v1 = 0
```

## Selector-output target

Freeze one explicit preobserver selector-output carrier:

```text
Q_out_v1 := span{o_+, o_-}
```

with ordered basis:

```text
o_+, o_-
```

## Exported selector output operator

Define the first actual selector-output map:

```text
O_sel_preLM_v1 : Q_sel_v1 -> Q_out_v1
```

by

```text
O_sel_preLM_v1(q_+) := o_+
O_sel_preLM_v1(q_-) := o_-
```

Equivalently, in the ordered bases `[q_+, q_-] -> [o_+, o_-]`:

```text
O_sel_preLM_v1 =
[[1, 0],
 [0, 1]]
```

## Source-side selector output

For the current source object:

```text
O_sel_preLM_v1(R_sel_preLM_v1(S_preLM_strict_core_source_object_v1))
  = r_plus_v1 o_+ + r_minus_v1 o_-
```

So the output-side response preserves:

```text
o_plus_v1  = r_plus_v1  > 0
o_minus_v1 = r_minus_v1 = 0
```

## Why this is an honest output step

1. it is derived only from `R_sel_preLM_v1`,
2. it stays preobserver,
3. it exports an actual selector-output map instead of only a reduction map,
4. it uses no imported `psi0`,
5. it uses no observer information,
6. it uses no external selector control,
7. it keeps `K_strict_gate` at operational-coefficient scope only.

## Bridge-ready output

Freeze:

```text
preobserver_output_state_v1
  := O_sel_preLM_v1(R_sel_preLM_v1(S_preLM_strict_core_source_object_v1))
```

This is only a future emergent-observer-side start state.

It does not export an actual emergent observer.

## Hard limits

`F91` does not claim:

- actual emergent observer construction,
- downstream completion beyond the output operator,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is:

```text
test whether O_sel_preLM_v1 is already an admissible strict-core preobserver
selector output operator
```
