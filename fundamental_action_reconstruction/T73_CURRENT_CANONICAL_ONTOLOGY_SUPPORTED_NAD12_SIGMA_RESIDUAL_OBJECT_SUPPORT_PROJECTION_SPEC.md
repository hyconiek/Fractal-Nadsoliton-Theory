# T73 Current Canonical-Ontology-Supported Nad12-Sigma Residual Object-Support Projection Spec

Status: `T73_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_RESIDUAL_OBJECT_SUPPORT_PROJECTION_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N337`, the nad12-sigma residual route already exports:

```text
one actual packaged residual object-support refinement candidate
```

The next honest question is narrower:

```text
can the current repo now honestly export one actual
residual object-support projection layer for this route,
using the newly exported route material above `N302`,
without pretending that any actual object support or
actual bridge/export-map object support is already present?
```

`T73` is intentionally weaker than:

```text
actual nad12-sigma object support
actual residual bridge/export-map object support
actual feeder support
actual theta export
actual pair population
actual loop break
actual E_orient
strict-core selector closure
ToE closure
```

## Scope

`T73` is scoped only to the following already exported lanes:

1. `N330`
   - actual object-support carrier/projection layer,
2. `N331`
   - actual packaged object-support witness candidate,
3. `N336`
   - actual packaged feeder-support candidate,
4. `N337`
   - actual packaged residual object-support refinement candidate,
5. `N299`
   - actual residual bridge-map target support,
6. `N302`
   - actual residual bridge/export-map object support still frozen below
     discharge,
7. sandbox `N18`
   - the theta/population loop still not actually broken.

## Residual object-support projection form

The strongest admissible projection form at this stage, if any, is only:

```text
residual_object_support_projection := {
  provider_route: "nad12_sigma_residual",
  object_support_carrier_projection: present,
  object_support_witness_candidate: present,
  feeder_support_candidate: present,
  residual_object_support_refinement_candidate: present,
  residual_target_support: present,
  status: "actual_projection_below_object_support"
}
```

with the exact reading:

1. the route may now be projected into the residual object-support frontier
   using newly accumulated route-local material,
2. that projection remains below actual nad12-sigma object support,
3. that projection remains below actual residual bridge/export-map object
   support,
4. no actual feeder support is exported,
5. no actual `theta_1`, `theta_2` are exported,
6. no actual pair population is exported,
7. no actual loop break is exported.

## Exact decision question

Can the current repo now honestly export anything stronger than:

```text
actual packaged residual object-support refinement candidate
```

namely:

```text
one actual residual object-support projection layer
for the nad12-sigma residual route,
using the already exported carrier/projection layer,
object-support witness candidate,
feeder-support candidate,
residual object-support refinement candidate,
and residual target-support layer,
while remaining explicitly below:
- actual object support,
- actual residual bridge/export-map object support,
- actual feeder support,
- actual theta export,
- actual pair population,
- actual loop break?
```

## Target

If the answer is positive at projection-only level, export:

```text
Xi_nad12_sigma_residual_object_support_projection_v1
```

with the intended meaning:

```text
actual residual object-support projection layer
for a future nad12-sigma object-support packet
and for a future residual bridge/export-map object-support packet

above refinement-candidate language
below actual object support
below actual residual bridge/export-map object support
```

## Hard limits

`T73` must not claim:

1. actual nad12-sigma object support,
2. actual residual bridge/export-map object support,
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
