# T120 Current Strict Sigma-Int Residual Bridge/Export-Map Object-Support Projection Spec

Status: `T120_CURRENT_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`N302` isolates the residual-datum / `sigma_int_candidate` third-provider route
below actual bridge/export-map object support.

`N384` now exports one corridor-protected **candidate** object-to-map support
projection artifact, i.e. a noncyclic, observer-free, pair-indexed mapping from
`sigma_int_input ∈ {+1,-1}` into a residual-datum target-slot *candidate*
population record (with admissible instantiations including
`sigma_int_candidate` and `sigma_int_strict_derived_v1`, without implying any
identification theorem between those two objects).

The next honest question is therefore narrower than closure:

```text
can the current repo export one actual route-level
residual bridge/export-map object-support projection layer
for the sigma-int third-provider route,
while remaining explicitly below:
- actual bridge/export-map object support,
- actual theta source,
- admissible S_sel_int,
- selector closure,
- ToE closure?
```

`T120` is intentionally weaker than any claim of discharge.

## Scope

`T120` is scoped only to the following already exported lanes:

1. `N299`
   - actual residual bridge-map target support is present,
2. `F311/N422`
   - the strict sigma-int lane exports an actual strict-core export-map object
     satisfying `T148` (residual `Z2` population only),
3. `N300/N301`
   - historical map-layer boundary/target (superseded/discharged by `F311/N422`),
4. `N302`
   - actual bridge/export-map object support remains frozen below discharge,
5. `R1`
   - strict-core residual-datum target-slot export object remains present,
6. `N384`
   - one corridor-protected object-to-map support projection **candidate**
     artifact exists (noncyclic + observer-free).

## Residual bridge/export-map object-support projection form

The strongest admissible projection form at this stage, if any, is only:

```text
residual_datum_sigma_int_bridge_export_map_object_support_projection := {
  provider_route: "sigma_int_third_provider_residual",
  residual_target_support: present,
  export_map_object_exported: present,
  export_map_nonexport_boundary: historical_superseded,
  export_map_object_target: historical_discharged,
  residual_target_slot_export: present,
  object_to_map_support_projection_candidate: present,
  status: "actual_projection_below_bridge_export_map_object_support"
}
```

with the exact reading:

1. the route may now be projected into the residual bridge/export-map
   object-support frontier using the corridor-protected projection candidate,
2. that projection remains below actual bridge/export-map object support,
3. no bridge/export map is exported,
4. no actual theta source is exported,
5. no admissible `S_sel_int` is exported,
6. no selector closure or `QW-2191` discharge is implied.

## Target

If the answer is positive at projection-only level, export:

```text
Xi_residual_datum_sigma_int_bridge_export_map_object_support_projection_v1
```

with intended meaning:

```text
actual residual bridge/export-map object-support projection layer
for the sigma-int third-provider route

above: object-target-only + template-carrier-only narrative
below: actual bridge/export-map object support
```

## Hard limits

`T120` must not claim:

1. discharge of `T2`,
2. actual residual bridge/export-map object support (the `N302` boundary is not
   discharged here),
3. actual `theta_1`, `theta_2`,
4. actual pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. `QW-2191` discharge,
8. ToE closure.
