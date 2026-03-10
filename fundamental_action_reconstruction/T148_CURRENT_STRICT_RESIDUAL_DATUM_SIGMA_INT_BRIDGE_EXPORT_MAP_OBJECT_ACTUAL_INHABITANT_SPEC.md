# T148 Current Strict Residual Datum Sigma-Int Bridge/Export-Map Object Actual Inhabitant Spec

Status: `T148_CURRENT_STRICT_RESIDUAL_DATUM_SIGMA_INT_BRIDGE_EXPORT_MAP_OBJECT_ACTUAL_INHABITANT_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

The strict-core residual-datum sigma-int bridge lane already exports:

1. a future-only target object naming the missing bridge/export-map object
   (`T36/F190/P281/N301`), and
2. a theorem-level map-layer nonexport boundary (`T35/F189/P280/N300`).

However, the lane still lacks a clear *discharge construction* spec for what
would count as exporting that missing bridge/export-map object **as an actual
strict-core object**, while keeping the strict discipline:

- no silent candidate→source upgrade,
- no silent gauge fixing,
- no silent overlay→selector-track promotion,
- no implied `QW-2191` discharge,
- no implied ToE closure.

`T148` provides that missing discharge construction spec in the same style as
`T145/T146`, but scoped to the sigma-int bridge/export-map object.

## Scope

`T148` is scoped only to the strict-core object:

```text
strict_core_equivalence_or_export_map : sigma_int_candidate -> residual orientation datum
```

It does not decide:

1. actual theta export,
2. actual pair population,
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.

## Output object

Export one **actual** strict-core bridge/export-map object:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_v1
```

intended to discharge the existing future-only target object:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_target_v1
```

(`T36/F190/N301`).

## Intended type-shape

`T148` does not mandate the internal implementation language of the map.
It only requires an explicit typed map-shape of the form:

```text
E_sigma_int_to_residual_datum_bridge_export_map_object_v1 :
  sigma_int_candidate -> residual_orientation_datum_target_slot
```

where `residual_orientation_datum_target_slot` is the strict-core slot already
materialized by `R1`.

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Upsilon_residual_datum_sigma_int_bridge_export_map_object_target_v1`
must at minimum provide:

1. **Actual map object export:** an exported strict-core object
   `Upsilon_residual_datum_sigma_int_bridge_export_map_object_v1` that carries
   an explicit mapping contract from `sigma_int_candidate` into the `R1`
   residual target slot.
2. **Strict source-upgrade discipline:** the export must not silently treat
   `sigma_int_candidate` as strict-derived. It must either:
   - discharge `T124/N389` by exporting a strict derivation/source-upgrade
     witness for sigma-int, or
   - keep sigma-int explicitly marked as candidate-only and therefore keep the
     map scoped strictly *below* strict-core bridge discharge (in which case
     the target is **not** discharged).
3. **Gauge-quotient discipline:** the export must not silently assume gauge
   safety. It must either:
   - discharge `T123/N388` by exporting a theorem-level gauge-quotient safety
     witness, or
   - keep gauge safety explicitly marked open (in which case the target is
     **not** discharged).
4. **Selector-track identification discipline:** the export must not silently
   promote overlay compatibility (`B7`) into selector-track identification. It
   must either:
   - discharge `T147/N414` by exporting a strict selector-track identification
     witness, or
   - keep selector-track identification explicitly marked overlay-only (in
     which case the target is **not** discharged).
5. **Noncyclic contract:** the map object must not take as input:
   - `theta_{1,2}`,
   - a populated basis-pair instance,
   (respects `N18`).
6. **Observer-free contract:** the map object must not use `K_obs`-indexed
   selection as a source of uniqueness.
7. **No implied selector closure:** the export must explicitly state:
   - `admissible_S_sel_int_present = false`,
   - `strict_core_selector_closure_present = false`,
   - `QW_2191_discharged = false`,
   unless those are separately discharged by new exported objects.

## Hard limits

`T148` must not claim:

1. discharge without explicitly addressing items (2)–(4) above,
2. discharge of `N300` without exporting an actual map object,
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.

