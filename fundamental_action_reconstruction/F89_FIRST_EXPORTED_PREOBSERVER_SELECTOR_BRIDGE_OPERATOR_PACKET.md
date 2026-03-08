# F89 First Exported Preobserver Selector Bridge Operator Packet

Status: `F89_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SELECTOR_BRIDGE_OPERATOR_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N196`, the next honest constructive move is no longer another negative
branch split.

It is one actual upstream selector-bridge operator built directly from the
exported preobserver orientation datum:

```text
E_orient_preLM_v1 -> B_sel_preLM_v1
```

without pretending that `R_sel`, `O_sel`, selector closure, or `QW-2191`
discharge already exist.

## Reused admissible orientation datum

Reuse:

```text
E_orient_preLM_v1 := (e_parallel, e_transverse)
```

on the topological-light plane `span{u_T, u_L}`, where:

```text
e_parallel   = (u_T + cos(phi) u_L) / sqrt(1 + cos(phi)^2)
e_transverse = (-cos(phi) u_T + u_L) / sqrt(1 + cos(phi)^2)
```

## Exported selector bridge operator

Freeze the first actual preobserver selector bridge operator as:

```text
B_sel_preLM_v1 := |e_parallel><e_parallel| - |e_transverse><e_transverse|
```

on `span{u_T, u_L}`, extended by zero on the matter slot `u_M`.

In the ordered basis `[u_T, u_L, u_M]`, let:

```text
c   := cos(phi)
N^2 := 1 + c^2
a   := (1 - c^2) / N^2
b   := 2 c / N^2
```

Then:

```text
B_sel_preLM_v1 =
[[ a,  b, 0],
 [ b, -a, 0],
 [ 0,  0, 0]]
```

## Immediate selector decomposition

Freeze the selector projectors:

```text
P_sel_plus_v1  := (I_TL + B_sel_preLM_v1) / 2 = |e_parallel><e_parallel|
P_sel_minus_v1 := (I_TL - B_sel_preLM_v1) / 2 = |e_transverse><e_transverse|
```

where `I_TL` is the identity on `span{u_T, u_L}`, extended by zero on `u_M`.

## Source-alignment witness

For

```text
S_preLM_strict_core_source_object_v1
  = u_T + cos(phi) u_L + (cos(phi)/4) u_M
```

its topological-light projection is:

```text
pi_TL(S_preLM_strict_core_source_object_v1) = u_T + cos(phi) u_L
                                           = sqrt(1 + cos(phi)^2) e_parallel
```

Therefore:

```text
P_sel_plus_v1  pi_TL(S_preLM_strict_core_source_object_v1)
  = pi_TL(S_preLM_strict_core_source_object_v1)

P_sel_minus_v1 pi_TL(S_preLM_strict_core_source_object_v1)
  = 0
```

and the signed selector response is strictly positive:

```text
<pi_TL(S), B_sel_preLM_v1 pi_TL(S)> = 1 + cos(phi)^2 > 0
```

## Why this is an honest bridge step

1. it is derived only from `E_orient_preLM_v1`,
2. it is upstream of observer,
3. it uses no imported `psi0`,
4. it uses no observer information,
5. it uses no external selector control,
6. it keeps `K_strict_gate` only at operational-coefficient scope,
7. it creates an actual signed selector decomposition instead of only a frame.

## Hard limits

`F89` does not claim:

- actual `R_sel`,
- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is:

```text
test whether B_sel_preLM_v1 is already an admissible strict-core preobserver
selector bridge operator
```
