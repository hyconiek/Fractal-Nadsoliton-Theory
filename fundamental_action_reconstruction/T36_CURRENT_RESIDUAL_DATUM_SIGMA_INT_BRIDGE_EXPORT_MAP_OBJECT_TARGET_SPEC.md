# T36 Current Residual Datum Sigma-Int Bridge Export Map Object Target Spec

Status: `T36_CURRENT_RESIDUAL_DATUM_SIGMA_INT_BRIDGE_EXPORT_MAP_OBJECT_TARGET_SPEC_NO_FALSE_PASS`
As of: `2026-03-11`

## Goal

Before `F311/N422`, the residual-datum / `sigma_int_candidate` route was
sharply frozen (via `N300/N301`) on the fact that:

```text
the actual bridge/export map was not yet exported in strict core
```

On the current repo state (`P388/P391`), the strict sigma-int lane now exports
an **actual** strict-core bridge/export-map object satisfying `T148`
(`F311/N422`).

Therefore `T36` is no longer a “current missing-object naming” spec.
It is kept as:

1. a sharp historical target-spec for the map-object shape, and
2. an admissible target-name record (now superseded by the actual export).

`T36` does **not** claim anything about post-`T148` object support, theta
export, selector closure, or ToE closure.

## Scope

`T36` is scoped only to the bridge/export-map object shape and its historical
target naming:

```text
strict_core_equivalence_or_export_map
  : sigma_int_candidate -> residual orientation datum
```

using only the current route-local material.

It does **not** decide:

1. actual bridge/export-map object support above the map object,
2. actual theta-source export,
3. actual component-2 support,
4. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
5. ToE closure,
6. all future bridge objects in general.

## Reused support

`T36` may reuse only already exported material:

1. `P4/P5`
   - the missing object is already sharply named as:
     `strict-core equivalence/export map sigma_int_candidate -> residual orientation datum`,
2. `T2`
   - one conditional bridge theorem spec already requires exactly that object,
3. `B4`
   - `sigma_int_candidate` exists as one candidate source object,
4. `C37`
   - the residual datum slot and candidate-fit already exist,
5. `R1`
   - the residual target-slot export object already exists,
6. `F311/N422`
   - an actual strict-core bridge/export-map object satisfying `T148` is now
     exported (discharging the earlier missing-object claim),
7. `C40/C41`
   - minimal field list and acceptance artifact schema already exist,
8. `C42/C43/C44/C45/C46`
   - dedicated carrier grammar, template content, admission, and minimal
     persisted carrier file already exist,
9. `F188/P279/N299`
   - the route already has actual support for the bridge-map layer,
10. `T35/F189/P280/N300`
   - the pre-`T148` map-layer nonexport boundary (historical).

## Exact decision question

On the current repo state, the repo exports strictly stronger material than:

```text
one explicit future-only target object is now named for the missing
strict-core bridge/export map,
```

because the map object is now actually exported (`F311/N422`).

So the only remaining role of `T36` is as a retained target-spec / target-name
record.

## Target

`T36` exports the historical future-only target name:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_target_v1
```

This target name is now superseded by the actual exported map object:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_v1
```

exported by `F311/N422`.

Historical intended meaning (pre-`T148`):

```text
the current repo now names one exact future-only target object for the missing
strict-core bridge/export map from sigma_int_candidate to the residual
orientation datum,
because the source object, target slot, candidate-fit, theorem-spec
requirements, and carrier grammar are all already sharply localized,
but the bridge/export map object itself is exported on the updated repo state
(`F311/N422`), so this “future-only target” description is historical and
superseded
```

## Hard limits

`T36` must not claim:

1. that bridge/export-map object export implies actual object support above the
   map object,
2. actual theta-source export,
3. actual component-2 support,
4. actual `theta_1`, `theta_2`,
5. actual populated basis-pair instance,
6. actual `E_orient`,
7. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
8. ToE closure.
