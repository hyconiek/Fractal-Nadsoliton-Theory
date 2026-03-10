# P388 Current Strict Residual Datum Sigma-Int Bridge/Export-Map Object Actual-Inhabitant Discharge Probe

Status: `P388_EXECUTED_CURRENT_STRICT_RESIDUAL_DATUM_SIGMA_INT_BRIDGE_EXPORT_MAP_OBJECT_ACTUAL_INHABITANT_DISCHARGE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Probe whether the current repo already exports an **actual** strict-core
bridge/export-map object discharging:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_target_v1 (T36/N301),
```

under the discharge acceptance tests of `T148`.

## Background

The lane already exports:

1. an explicit map-layer nonexport boundary (`T35/F189/P280/N300`),
2. an explicit future-only target naming the missing bridge/export-map object
   (`T36/F190/P281/N301`),
3. stronger post-witness support on the provider-object witness track via an
   actual object-support carrier (`T146/F301/N413`).

`P388` checks whether any of that material already constitutes an **actual**
export-map object, or whether the bridge/export-map object still remains
absent exactly as stated by `N300`.

## Probe table (T148 acceptance tests)

| Acceptance test (T148) | Verdict | Evidence |
|---|---|---|
| actual map object `Upsilon_*_export_map_object_v1` exported | NO | `N300` boundary remains in force; no exported strict-core map object is present |
| typed map-shape `sigma_int_candidate -> residual_orientation_datum_target_slot` exported | NO | no exported strict-core object provides that contract |
| sigma-int strict derivation/source upgrade discharged | NO | only future-only target `T124/N389` exists |
| sigma-int gauge-quotient safety discharged | NO | only future-only target `T123/N388` exists |
| selector-track identification beyond overlay-only discharged | NO | only future-only target `T147/N414` exists |
| noncyclic contract (no theta/populated inputs) satisfied | N/A | map object absent |
| observer-free contract satisfied | N/A | map object absent |
| no implied selector closure / `QW-2191` discipline | N/A | map object absent |

## Exact verdict

The strict-core bridge/export-map object remains absent:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_target_v1: NOT DISCHARGED
N300 map-layer nonexport boundary: STILL IN FORCE
```

`N413` strengthens post-witness object-support on the provider-object witness
track but does **not** export an equivalence/export-map object.

## Consequence (next honest step)

The next honest move on this sublane is not to claim “bridge discharge”.

It is to discharge at least one of the strict prerequisites explicitly listed
by `T148` (on a declared strict domain), and only then attempt an actual map
object export:

1. `T124/N389` sigma-int strict derivation/source upgrade,
2. `T123/N388` sigma-int gauge-quotient safety,
3. `T147/N414` selector-track identification beyond overlay-only.

