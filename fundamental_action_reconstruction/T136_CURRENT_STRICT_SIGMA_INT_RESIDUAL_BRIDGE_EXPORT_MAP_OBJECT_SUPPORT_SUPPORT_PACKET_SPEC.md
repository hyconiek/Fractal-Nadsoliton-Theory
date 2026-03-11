# T136 Current Strict Sigma-Int Residual Bridge/Export-Map Object-Support Support-Packet Spec

Status: `T136_CURRENT_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_SUPPORT_PACKET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `N387`, the sigma-int residual third-provider route exports:

1. an object-support projection layer into the object-support frontier (`N385`),
2. an object-support witness candidate and witness (`N386`, `N387`),
3. a sharply named future-only target object for the next missing layer
   above the witness (`N395` / `T130`),

while remaining explicitly below:

- actual bridge/export-map object support (`N302`),
- strict-core export-map object support above the exported map object (still absent),
- discharge of `T2`,
- strict-core theta export / pair population (`N18` loop not broken),
- admissible `S_sel_int` / selector closure / `QW-2191` discharge,
- ToE closure.

The next honest move is **not** to relabel the route as discharged.

It is narrower:

```text
package one actual support packet for the next missing post-witness layer
(actual object support),
so that the post-N387 frontier is not only "missing",
but also carries one explicit support-packet object below N302.
```

`T136` specifies that support-packet form.

## Scope

`T136` is scoped only to:

1. `N299` (bridge-map target support),
2. `F311/N422` (export-map object exported; residual `Z2` population only),
3. `N300/N301` (historical map-layer boundary/target; superseded/discharged),
4. `N385` (object-support projection layer),
5. `N387` (object-support witness layer),
6. `N395` (future-only actual object-support target naming),
7. `N302` (no actual object support on the current repo state),
8. convergence-side coherence support (optional strengthening only):
   - `N401` (joint coherence witness candidate),
   - `N402` (pullback support carrier candidate),
   both explicitly below `N302`.

## Support-packet form (admissible at this stage)

The strongest admissible support-packet form at this stage, if any, is only:

```text
sigma_int_residual_bridge_export_map_object_support_support_packet := {
  provider_route: "sigma_int_third_provider_residual",

  bridge_map_target_support: present,
  export_map_object_exported: present,
  export_map_nonexport_boundary: historical_superseded,
  export_map_object_target: historical_discharged,

  object_support_projection_layer: present,
  object_support_witness: present,
  actual_object_support_target_named: future_only_present,

  convergence_side_joint_coherence_witness_candidate: optional_present,
  convergence_side_pullback_support_carrier_candidate: optional_present,

  status: "actual_support_packet_below_bridge_export_map_object_support"
}
```

with the exact reading:

1. the route may now carry one actual support packet for its next missing
   post-witness layer (actual object support),
2. that support packet remains below actual object support (`N302`),
3. no strict-core theta export / pair population is implied,
4. no selector closure is implied.

## Target

If the answer is positive at support-packet level, export:

```text
Nu_residual_datum_sigma_int_bridge_export_map_object_support_support_packet_v1
```

with intended meaning:

```text
actual support packet for the next missing post-witness layer
on the sigma-int residual third-provider route

above: witness-only language (N387)
below: actual bridge/export-map object support (N302)
```

## Hard limits

`T136` must not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. discharge of `T2`,
4. actual theta export / pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.
