# T138 Current Strict ToE Closure Provider-Object Carrier Residual Bridge/Export-Map Object-Support Witness Spec

Status: `T138_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_WITNESS_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `T137`, the strict provider-object carrier residual branch may already
export:

```text
one actual packaged bridge/export-map object-support witness candidate
```

The next honest question is narrower:

```text
can the current repo now honestly export
one actual bridge/export-map object-support witness
on the same provider-object carrier route,
without pretending that actual bridge/export-map object support
or bridge/export map export
is already present?
```

`T138` is intentionally weaker than:

```text
actual bridge/export-map object support
bridge/export map export
strict-core theta export
pair population
admissible S_sel_int
selector closure
ToE closure
```

## Scope

`T138` is scoped only to the already exported provider-object carrier residual
bridge lane:

1. `N299` (bridge-map target support),
2. `F311/N422` (actual export-map object satisfying `T148`; supersedes `N300`
   and discharges `N301`),
3. `N300/N301` (historical map-layer boundary/target; superseded/discharged),
4. `N397` (weld discharge at declared interface level),
5. `N399` (corridor-protected provider-carrier projection candidate),
6. `N398` (actual projection layer into object-support frontier),
7. `T137` (witness-candidate spec only).

## Witness form (admissible at this stage)

The strongest admissible witness form at this stage, if any, is only:

```text
provider_object_carrier_residual_bridge_export_map_object_support_witness := {
  provider_route: "provider_object_carrier_orbit_quotient",
  weld_discharge_present: true,
  bridge_map_target_support: present,
  export_map_nonexport_boundary: historical_only,
  export_map_object_target: discharged,
  export_map_object_exported: true,
  object_to_map_support_projection_candidate: present,
  object_support_projection_layer: present,
  object_support_witness_candidate: present,
  status: "actual_witness_below_bridge_export_map_object_support"
}
```

with the exact reading:

1. the route may now carry one actual witness for the missing bridge/export-map
   object-support layer,
2. that witness remains below any actual bridge/export-map object support,
3. no bridge/export map is exported.

## Target

If the answer is positive at witness level, export:

```text
Kappa_residual_datum_provider_object_carrier_bridge_export_map_object_support_witness_v1
```

with the intended meaning:

```text
actual bridge/export-map object-support witness

above witness-candidate-only language
below actual bridge/export-map object support
below bridge/export map export
```

## Hard limits

`T138` must not claim:

1. discharge of `T2`,
2. actual bridge/export-map object support (`N395` remains future-only),
3. any new bridge/export-map object export beyond `F311/N422`,
4. discharge of `N395`,
5. actual `theta_1`, `theta_2`,
6. actual pair population,
7. admissible `S_sel_int`,
8. strict-core selector closure,
9. `QW-2191` discharge,
10. ToE closure.
