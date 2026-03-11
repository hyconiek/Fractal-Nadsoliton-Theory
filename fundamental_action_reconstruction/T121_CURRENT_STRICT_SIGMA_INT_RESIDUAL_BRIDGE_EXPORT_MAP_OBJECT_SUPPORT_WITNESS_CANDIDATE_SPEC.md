# T121 Current Strict Sigma-Int Residual Bridge/Export-Map Object-Support Witness Candidate Spec

Status: `T121_CURRENT_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_WITNESS_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`N385` exports one actual projection layer of the residual-datum /
`sigma_int_candidate` third-provider route into the bridge/export-map
object-support frontier, while remaining explicitly below actual object support.

The next honest positive move, if any, is still not discharge.

It is to package the route as one explicit **object-support witness candidate**
for that missing bridge/export-map object-support layer, without claiming:

1. actual bridge/export-map object support,
2. export of the bridge/export map,
3. strict-core actual `theta_1,theta_2`,
4. admissible `S_sel_int`,
5. selector closure or `QW-2191` discharge,
6. ToE closure.

## Scope

`T121` is scoped only to the already exported sigma-int residual bridge lane:

1. `N299` (bridge-map target support),
2. `F311/N422` (export-map object exported; residual `Z2` population only),
3. `N300/N301` (historical map-layer boundary/target; superseded/discharged),
4. `N384` (corridor-protected projection candidate artifact),
5. `N385` (actual projection layer into the object-support frontier),
6. `R1` (residual-datum target-slot export packet).

## Witness candidate form (admissible at this stage)

The strongest admissible witness-candidate form at this stage is only:

```text
residual_datum_sigma_int_bridge_export_map_object_support_witness_candidate := {
  provider_route: "sigma_int_third_provider_residual",
  residual_target_slot_export: present,
  bridge_map_target_support: present,
  export_map_object_exported: present,
  export_map_nonexport_boundary: historical_superseded,
  export_map_object_target: historical_discharged,
  object_to_map_support_projection_candidate: present,
  residual_object_support_projection_layer: present,
  status: "actual_witness_candidate_below_bridge_export_map_object_support"
}
```

with the exact reading:

1. the route is now strong enough to be packaged as a candidate witness for the
   missing bridge/export-map object-support layer,
2. but this does not export any actual bridge/export-map object support.

## Target

If the packaging is admissible, export:

```text
Omega_residual_datum_sigma_int_bridge_export_map_object_support_witness_candidate_v1
```

## Hard limits

`T121` must not claim:

1. discharge of `T2`,
2. actual bridge/export-map object support,
3. actual `theta_1`, `theta_2`,
4. actual pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. `QW-2191` discharge,
8. ToE closure.
