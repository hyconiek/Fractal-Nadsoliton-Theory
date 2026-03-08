# F90 First Exported Preobserver Selector Reduction Operator Packet

Status: `F90_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SELECTOR_REDUCTION_OPERATOR_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N197`, the next honest constructive move is:

```text
B_sel_preLM_v1 -> R_sel_preLM_v1
```

That is, export one actual preobserver selector-reduction map from the already
constructed selector bridge operator, without pretending that `O_sel`,
downstream completion, selector closure, or `QW-2191` discharge already exist.

## Reused selector bridge input

Reuse:

```text
B_sel_preLM_v1 := |e_parallel><e_parallel| - |e_transverse><e_transverse|
```

with selector projectors:

```text
P_sel_plus_v1  := |e_parallel><e_parallel|
P_sel_minus_v1 := |e_transverse><e_transverse|
```

on the typed preobserver carrier:

```text
V_topo ⊕ L_int ⊕ M_int
```

with basis `[u_T, u_L, u_M]`.

## Selector-sector reduction target

Freeze one explicit selector-sector output carrier:

```text
Q_sel_v1 := span{q_+, q_-}
```

with ordered basis:

```text
q_+, q_-
```

## Exported selector reduction map

Define the first actual selector reduction map:

```text
R_sel_preLM_v1 : V_topo ⊕ L_int ⊕ M_int -> Q_sel_v1
```

by

```text
R_sel_preLM_v1(x)
  := (<e_parallel, x>) q_+ + (<e_transverse, x>) q_-
```

Equivalently, in the ordered bases `[u_T, u_L, u_M] -> [q_+, q_-]`:

```text
R_sel_preLM_v1 =
[[ e_parallel^T ],
 [ e_transverse^T ]]
```

that is:

```text
R_sel_preLM_v1 =
[[  0.71177976495,  0.70240270942, 0.0 ],
 [ -0.70240270942,  0.71177976495, 0.0 ]]
```

## Source-side selector response

For

```text
S_preLM_strict_core_source_object_v1
  = u_T + cos(phi) u_L + (cos(phi)/4) u_M
```

its selector-sector reduction is:

```text
R_sel_preLM_v1(S_preLM_strict_core_source_object_v1)
  = sqrt(1 + cos(phi)^2) q_+ + 0 q_-
```

So the source-side reduction has:

```text
r_plus_v1  = sqrt(1 + cos(phi)^2) > 0
r_minus_v1 = 0
```

## Why this is an honest reduction step

1. it is derived only from `E_orient_preLM_v1 / B_sel_preLM_v1`,
2. it stays preobserver,
3. it exports an actual selector-sector map instead of only a bridge operator,
4. it uses no imported `psi0`,
5. it uses no observer information,
6. it uses no external selector control,
7. it keeps `K_strict_gate` at operational-coefficient scope only.

## Bridge-ready output

Freeze:

```text
O_sel_start_state_v1 := R_sel_preLM_v1(S_preLM_strict_core_source_object_v1)
```

This is only a future observer-side start state.

It does not export `O_sel`.

## Hard limits

`F90` does not claim:

- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is:

```text
test whether R_sel_preLM_v1 is already an admissible strict-core preobserver
selector reduction operator
```
