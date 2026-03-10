# T78 Current Canonical-Ontology-Supported Nad12-Sigma Residual Bridge/Export-Map Object-Support Witness Spec

Status: `T78_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_WITNESS_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N342`, the nad12-sigma residual route already exports:

```text
one actual support packet
for the next missing residual bridge/export-map object-support layer
```

The next honest question is narrower:

```text
can the current repo now honestly export
one actual residual bridge/export-map object-support witness
on the same route,
without pretending that actual residual bridge/export-map object support
or actual nad12-sigma object support
is already present?
```

`T78` is intentionally weaker than:

```text
actual residual bridge/export-map object support
actual nad12-sigma object support
actual feeder support
actual theta export
actual pair population
actual loop break
actual E_orient
strict-core selector closure
ToE closure
```

## Scope

`T78` is scoped only to the following already exported lanes:

1. `N299`
   - actual residual bridge-map target support,
2. `N302`
   - actual residual bridge/export-map object support still frozen below
     discharge,
3. `N338`
   - actual residual object-support projection,
4. `N339`
   - actual nad12-sigma residual object-support witness,
5. `N340`
   - actual nad12-sigma residual object-support support packet,
6. `N341`
   - actual residual bridge/export-map support witness,
7. `N342`
   - actual residual bridge/export-map object-support support packet,
8. sandbox `N18`
   - theta/population loop still not actually broken.

## Residual object-support witness form

The strongest admissible witness form at this stage, if any, is only:

```text
nad12_sigma_residual_bridge_export_map_object_support_witness := {
  provider_route: "nad12_sigma_residual",
  residual_target_support: present,
  residual_support_witness: present,
  residual_object_support_projection: present,
  residual_object_support_witness: present,
  residual_object_support_packet: present,
  residual_object_support_support_packet: present,
  status: "actual_witness_below_residual_object_support"
}
```

with the exact reading:

1. the route may now carry one actual witness for the missing residual
   bridge/export-map object-support layer,
2. that witness remains below actual residual bridge/export-map object
   support,
3. that witness remains below actual nad12-sigma object support,
4. no actual feeder support is exported,
5. no actual `theta_1`, `theta_2` are exported,
6. no actual pair population is exported,
7. no actual loop break is exported.

## Exact decision question

Can the current repo now honestly export anything stronger than:

```text
actual support packet
for the next missing residual bridge/export-map object-support layer
```

namely:

```text
one actual residual bridge/export-map object-support witness
on the same nad12-sigma residual route,
using the already exported residual target-support layer,
residual bridge/export-map support witness,
residual object-support projection,
residual object-support witness,
object-support support packet,
and residual bridge/export-map object-support support packet,
while remaining explicitly below:
- actual residual bridge/export-map object support,
- actual nad12-sigma object support,
- actual feeder support,
- actual theta export,
- actual pair population,
- actual loop break?
```

## Target

If the answer is positive at witness level, export:

```text
Tau_nad12_sigma_residual_bridge_export_map_object_support_witness_v1
```

with the intended meaning:

```text
actual residual bridge/export-map object-support witness

above support-packet-only language
below actual residual bridge/export-map object support
below actual nad12-sigma object support
```

## Hard limits

`T78` must not claim:

1. actual residual bridge/export-map object support,
2. actual nad12-sigma object support,
3. actual sigma asymmetry law,
4. actual `f(sigma_int)` weighting law,
5. actual `lambda_1`, `lambda_2`,
6. actual `u_1`, `u_2`,
7. actual `theta_1`, `theta_2`,
8. actual populated basis-pair instance,
9. actual loop break,
10. actual `E_orient`,
11. admissible `S_sel_int`,
12. strict-core selector closure,
13. `QW-2191` discharge,
14. ToE closure.
