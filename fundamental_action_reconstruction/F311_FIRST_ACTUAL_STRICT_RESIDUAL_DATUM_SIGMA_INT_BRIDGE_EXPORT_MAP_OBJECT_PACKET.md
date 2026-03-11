# F311 First Actual Strict Residual Datum Sigma-Int Bridge/Export-Map Object Packet

Status: `F311_EXECUTED_FIRST_ACTUAL_STRICT_RESIDUAL_DATUM_SIGMA_INT_BRIDGE_EXPORT_MAP_OBJECT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `P388`, the strict-core residual-datum sigma-int bridge lane was blocked
by one missing object-layer item:

```text
an actual strict-core bridge/export-map object discharging N301 (T148).
```

The upstream strict prerequisites demanded by `T148` are now exported:

1. strict sigma-int source upgrade (`F307/N418`),
2. theorem-level gauge-quotient safety (`F308/N419`),
3. selector-track identification beyond overlay-only (`F310/N421`),

all while explicitly keeping `QW-2191` open.

This packet executes the next honest move:

```text
export one actual strict-core bridge/export-map object
Upsilon_residual_datum_sigma_int_bridge_export_map_object_v1
with explicit typed map-shape into the R1 residual target slot,
without theta inputs and without implied selector closure.
```

## Inputs reused (strict-admissible)

1. `R1`
   - strict-core residual orientation datum target-slot export.
2. `T148`
   - discharge acceptance tests for an actual map object.
3. `F307/N418`
   - strict sigma-int source upgrade exported.
4. `F308/N419`
   - theorem-level gauge-quotient safety exported.
5. `F310/N421`
   - selector-track identification witness exported (keeps `QW-2191` open).

## Exported map object

Export one strict-core bridge/export-map object:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_v1.
```

### Typed map-shape

Export an explicit typed map-shape:

```text
E_sigma_int_to_residual_datum_bridge_export_map_object_v1 :
  sigma_int_strict_derived_v1 -> residual_orientation_datum_target_slot
```

where:

1. `sigma_int_strict_derived_v1` is the strict-core sigma-int source-upgrade
   datum exported in `F307`,
2. `residual_orientation_datum_target_slot` is the strict-core target slot
   exported in `R1`.

### Map meaning (residual Z2 population only; no theta export)

The exported map object populates **only** the residual `Z2` sign convention
layer of the residual orientation datum target slot, by the rule:

```text
residual_orientation_sign_convention := sigma_int_strict_derived_v1 ∈ {+1,-1}.
```

Implementation note (purely formal, no theta inputs):

Because the `R1` slot is represented by the basis formulas
`u_1(theta_1), u_2(theta_2)`, the residual `Z2` flip can be represented by the
formal basis sign flip:

```text
u_1  ->  (sigma_int_strict_derived_v1) * u_1,
u_2  ->  u_2,
```

equivalently (when desired) by the formal shift:

```text
theta_1 -> theta_1 + π    when sigma_int_strict_derived_v1 = -1,
```

without taking `theta_1, theta_2` as inputs and without claiming any strict
theta-source export.

## Persisted artifact

`fundamental_action_reconstruction/generated/upsilon_residual_datum_sigma_int_bridge_export_map_object_v1.json`

## Status discipline

This packet does **not** claim:

1. export of strict-core theta sources (`theta_{1,2}`) or population of the
   residual target slot,
2. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
3. ToE closure.

