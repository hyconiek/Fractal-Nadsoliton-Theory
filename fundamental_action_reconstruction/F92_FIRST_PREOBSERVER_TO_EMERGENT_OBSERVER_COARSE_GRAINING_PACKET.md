# F92 First Preobserver To Emergent Observer Coarse Graining Packet

Status: `F92_EXECUTED_FIRST_PREOBSERVER_TO_EMERGENT_OBSERVER_COARSE_GRAINING_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N199`, the next honest constructive move is:

```text
O_sel_preLM_v1 -> C_obs_limit_preLM_v1
```

That is, export one actual downstream coarse-graining map from the already
constructed strict-core preobserver selector-output state into a macroscopic
observer-limit carrier, without pretending that an actual emergent observer,
selector closure, or `QW-2191` discharge already exist.

## Reused selector-output input

Reuse:

```text
O_sel_preLM_v1 : Q_sel_v1 -> Q_out_v1
```

with output basis:

```text
Q_out_v1 := span{o_+, o_-}
```

and current source-side output response:

```text
O_sel_preLM_v1(R_sel_preLM_v1(S_preLM_strict_core_source_object_v1))
  = o_plus_v1 o_+ + o_minus_v1 o_-
```

where:

```text
o_plus_v1  > 0
o_minus_v1 = 0
```

## Observer-limit coarse-graining target

Freeze one explicit downstream coarse-graining carrier:

```text
Y_obs_limit_v1 := span{y_bias, y_total}
```

with ordered basis:

```text
y_bias, y_total
```

Interpretation:

1. `y_bias` records selector asymmetry at macroscopic readout level,
2. `y_total` records total coarse-grained signal intensity,
3. both remain preobserver-to-observer-limit quantities,
4. neither quantity is yet an actual emergent observer state.

## Exported coarse-graining operator

Define the first explicit downstream coarse-graining map:

```text
C_obs_limit_preLM_v1 : Q_out_v1 -> Y_obs_limit_v1
```

by

```text
C_obs_limit_preLM_v1(o_+) := (1/2) y_bias + (1/2) y_total
C_obs_limit_preLM_v1(o_-) := (-1/2) y_bias + (1/2) y_total
```

Equivalently, in ordered bases `[o_+, o_-] -> [y_bias, y_total]`:

```text
C_obs_limit_preLM_v1 =
[[ 1/2, -1/2],
 [ 1/2,  1/2]]
```

## Source-side coarse-grained response

For the current source object:

```text
C_obs_limit_preLM_v1(
  O_sel_preLM_v1(R_sel_preLM_v1(S_preLM_strict_core_source_object_v1))
)
  = y_bias_v1 y_bias + y_total_v1 y_total
```

with:

```text
y_bias_v1  = (o_plus_v1 - o_minus_v1)/2
y_total_v1 = (o_plus_v1 + o_minus_v1)/2
```

So the current source produces:

```text
y_bias_v1  > 0
y_total_v1 > 0
```

because `o_plus_v1 > 0` and `o_minus_v1 = 0`.

## Why this is an honest downstream move

1. it is derived only from the already exported selector-output state,
2. it remains downstream of `nadsoliton -> light -> matter`,
3. it does not promote observer information into the primary selector source,
4. it provides a macroscopic observer-limit readout carrier without claiming
   an actual observer,
5. it uses no imported `psi0`,
6. it uses no external selector control,
7. it keeps `K_strict_gate` at operational-control scope only.

## Hard limits

`F92` does not claim:

- actual emergent observer construction,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is:

```text
test whether C_obs_limit_preLM_v1 is already an admissible strict-core
preobserver-to-emergent-observer coarse-graining operator
```
