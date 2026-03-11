# F301 First Actual Strict Residual Datum Bridge/Export-Map Object-Support Carrier Packet

Status: `F301_EXECUTED_FIRST_ACTUAL_STRICT_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_CARRIER_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Execute `T146` if admissible: export one actual post-witness object-support
carrier above the provider-object witness layer (`N405`), discharging the
carrier target from `T140/N407`, while remaining explicitly below **actual**
bridge/export-map object support above the map object (`N395`) and below
selector closure.

## Inputs reused

1. `T140/N407`
   - post-witness carrier target exists (future-only),
2. `N405`
   - provider-object witness exists,
3. `N410`
   - contraction parameters `(a,b)` are strict source-side outputs of `tau_src_candidate_v1`,
4. `N412`
   - actual provider-object carrier layer exists,
5. `T141/F296`
   - finite nad12-depth preobject carrier construction exists (candidate form),
6. `T146/P386`
   - carrier inhabitant spec and probe verdict.

## Packet result

`F301` exports one actual post-witness object-support carrier:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_v1
```

with structured content:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_v1 :=
{
  upgrades_witness:
    Kappa_residual_datum_provider_object_carrier_bridge_export_map_object_support_witness_v1,
  source_domain: tau_src_candidate_v1,
  contraction_parameter_source_map: A_strict_provider_object_contraction_parameter_source_map_v1,
  induced_parameters: (a,b) := (cos(phi), cos(phi)),
  carrier_space: ℓ^2(ℤ),
  gauge_symmetry_group: U(1),
  pair_index_set: {pair1, pair2},
  nad12_depth: K := 11,
  carrier_vectors: (ψ_support_pair1, ψ_support_pair2),
  bridge_facing_projection_layer:
    Xi_residual_datum_provider_object_carrier_bridge_export_map_object_support_projection_v1,
  export_map_nonexport_boundary: historical_only (N300; superseded by F311/N422),
  export_map_object_exported: true (F311/N422; discharges N301/T148),
  selector_neutral: true,
  ToE_closed: false
}.
```

## Persisted carrier artifact instance

`F301` also records one persisted carrier artifact instance:

```text
fundamental_action_reconstruction/generated/residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json
```

This file is a carrier artifact only. It is not an export-map object and does
not export any actual bridge/export-map object support above the map object
(`N395` remains future-only).

## Hard limits

`F301` must not claim:

1. any actual bridge/export-map object support above the map object (`N395`
   remains future-only),
2. discharge of sigma-int strict prerequisites (`T123/T124`),
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.
