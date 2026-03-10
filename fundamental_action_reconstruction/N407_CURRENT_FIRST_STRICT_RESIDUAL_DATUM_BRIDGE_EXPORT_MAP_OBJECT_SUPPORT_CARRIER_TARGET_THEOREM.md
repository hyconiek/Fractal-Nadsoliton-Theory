# N407 Current First Strict Residual Datum Bridge/Export-Map Actual Object-Support Carrier Target Theorem

Status: `N407_DISCHARGED_CURRENT_FIRST_STRICT_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_CARRIER_TARGET_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Package theorem-level the strongest honest current statement about the
post-witness **object-support carrier** layer on the strict residual-datum
bridge/export-map lane after the first route-local discharge instance
(`N413`), while staying explicitly below export-map export (`N300`), selector
closure, and ToE closure.

## Theorem-level conclusion

From `T140/P381/F295`, the current repo exports one future-only target object:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_target_v1
```

with the following exact meaning:

1. `N387` remains correct:
   - the sigma-int residual branch reaches one object-support witness layer,
2. `N405` remains correct:
   - the provider-object carrier residual branch reaches one object-support
     witness layer,
3. `N413` now exports one actual post-witness carrier:
   - `Omicron_residual_datum_bridge_export_map_object_support_carrier_v1`,
   - discharging the target on the provider-object witness track (`N405`),
   - without using `sigma_int_candidate` as a strict-core source datum,
   - while keeping `N300` (export-map nonexport) in force,
4. therefore the repo now has both:
   - sharp target naming (`T140`), and
   - one route-local discharge instance above witness (`N413`),
5. no bridge/export-map export, theta export, selector closure, or ToE closure
   claim is implied.

## What N407 proves

`N407` proves only this narrower statement:

1. the repo names one explicit target object for the post-witness object-support
   carrier layer (`T140`),
2. the repo exports one explicit discharge instance above the provider-object
   witness track (`N413`),
3. no bridge/export-map export, theta export, selector closure, or ToE closure
   claim is implied.

## What N407 does not prove

`N407` does not prove:

1. that the post-witness carrier discharges the sigma-int witness track
   (`N387`) without discharging `T123/T124`,
2. discharge of `N302` on the sigma-int third-provider route,
3. bridge/export-map export,
4. theta export / pair population,
5. admissible `S_sel_int`,
6. selector closure or `QW-2191` discharge,
7. ToE closure.

## Consequence (next honest step)

After `N407` (and `N413`), the next honest move is no longer “add a post-witness
carrier”.
It must be either:

1. break the export-map nonexport boundary by exporting one actual bridge/export
   map object (`N300`), or
2. introduce a genuinely new blocker-cut changing the map-layer diagnosis,

both noncyclically and observer-free, and without selector false-pass.
