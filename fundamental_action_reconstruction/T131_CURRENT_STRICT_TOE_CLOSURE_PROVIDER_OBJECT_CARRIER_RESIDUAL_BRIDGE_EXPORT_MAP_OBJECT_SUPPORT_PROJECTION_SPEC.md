# T131 Current Strict ToE Closure Provider-Object Carrier Residual Bridge/Export-Map Object-Support Projection Spec

Status: `T131_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After:

1. `T126/N391` (orbit-quotient provider-object carrier candidate exists),
2. `T127/N392` (bridge-facing projection candidate exists),
3. `N397` (welding target `T128` discharged at declared interface level),

the provider-object lane still does **not** export:

1. actual bridge/export-map object support (`N302` remains in force),
2. any new bridge/export-map object export beyond `F311/N422`,
3. actual theta export / pair population (`N18` remains in force),
4. admissible `S_sel_int` / selector closure (`QW-2191` remains in force).

The next honest move is therefore narrower than closure:

```text
can the current repo export one actual route-level
object-support projection layer for the provider-object carrier lane
into the residual bridge/export-map object-support frontier,
while remaining explicitly below N302?
```

This is the provider-object analogue of the sigma-int projection layer
packaged in `T120/N385`, but **without** claiming any sigma-int strict-source
upgrade.

## Scope

`T131` is scoped only to already exported material:

1. `T127/N392`
   - provider-object carrier → residual frontier projection candidate exists,
2. `T118/N384/N385`
   - downstream residual-frontier record class name is already stabilized on
     the sigma-int projection lane,
3. `R1`
   - residual-datum target-slot export language exists,
4. `N302`
   - object-support incompatibility boundary remains in force.

## Residual-frontier projection form (provider-object lane)

The strongest admissible provider-object projection form at this stage, if any,
is only:

```text
residual_datum_provider_object_carrier_bridge_export_map_object_support_projection := {
  provider_route: "provider_object_carrier_orbit_quotient",
  weld_discharge_present: true,
  object_to_map_support_projection_candidate: present,
  residual_frontier_record_class: "residual_datum_bridge_export_map_object_support_projection_candidate_instance",
  status: "actual_projection_layer_below_bridge_export_map_object_support"
}
```

## Target

If the answer is positive at projection-only level, export:

```text
Xi_residual_datum_provider_object_carrier_bridge_export_map_object_support_projection_v1
```

with intended meaning:

```text
actual residual bridge/export-map object-support projection layer
for the provider-object carrier lane

above: projection-candidate-only language (T127)
below: actual bridge/export-map object support (N302)
```

## Hard limits

`T131` must not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. any new export-map object export beyond `F311/N422`,
4. actual theta export / pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure / `QW-2191` discharge,
7. ToE closure.
