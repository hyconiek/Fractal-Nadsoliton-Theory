# F300 First Actual Strict ToE Closure Provider-Object Carrier Layer Packet

Status: `F300_EXECUTED_FIRST_ACTUAL_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_LAYER_PACKET_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Execute `T145` if admissible: export one actual strict provider-object carrier
layer inhabitant, discharging the carrier-layer target from `T125/N390`,
while staying explicitly below `N302` and below realization/selector closure.

## Inputs reused

1. `T125/N390`
   - provider-object carrier-layer target exists (future-only),
2. `T126/F279/N391`
   - orbit-quotient carrier form exists (candidate form),
3. `N410`
   - contraction parameters `(a,b)` are strict source-side outputs of `tau_src_candidate_v1`,
4. `N398`
   - bridge-facing residual-frontier projection layer exists (below `N302`),
5. `T145/P385`
   - discharge construction and probe verdict.

## Packet result

`F300` exports one actual carrier-layer inhabitant:

```text
Epsilon_strict_provider_object_carrier_layer_v1
```

with structured content:

```text
Epsilon_strict_provider_object_carrier_layer_v1 :=
{
  carrier_form: "orbit_quotient_on_l2Z",
  source_domain: tau_src_candidate_v1,
  contraction_parameter_source_map: A_strict_provider_object_contraction_parameter_source_map_v1,
  induced_parameters: (a,b) := (cos(phi), cos(phi)),
  pair_index_set: {pair1, pair2},
  gauge_symmetry_group: U(1),
  seed_state: ψ_src := δ_0,
  output_carrier: ([ψ_src]_pair1, [ψ_src]_pair2),
  bridge_facing_projection_layer: Xi_residual_datum_provider_object_carrier_bridge_export_map_object_support_projection_v1,
  n302_boundary: in_force,
  admissible_S_sel_int_present: false,
  strict_core_selector_closure_present: false,
  QW_2191_discharged: false,
  ToE_closed: false
}.
```

## Persisted carrier-layer artifact instance

`F300` also records one persisted carrier-layer artifact instance:

```text
fundamental_action_reconstruction/generated/provider_object_carrier_layer_actual_inhabitant_instance.json
```

This file is a carrier artifact only. It is not an export-map object and does
not discharge `N302`.

## Hard limits

`F300` must not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. any export-map object export,
4. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
5. provider-object realization,
6. ToE closure.

