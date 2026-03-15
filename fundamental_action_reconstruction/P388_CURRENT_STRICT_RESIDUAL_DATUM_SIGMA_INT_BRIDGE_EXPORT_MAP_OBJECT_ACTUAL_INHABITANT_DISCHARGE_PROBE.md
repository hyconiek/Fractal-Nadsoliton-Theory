# P388 Current Strict Residual Datum Sigma-Int Bridge/Export-Map Object Actual-Inhabitant Discharge Probe

Status: `P388_EXECUTED_CURRENT_STRICT_RESIDUAL_DATUM_SIGMA_INT_BRIDGE_EXPORT_MAP_OBJECT_ACTUAL_INHABITANT_DISCHARGE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Probe whether the current repo already exports an **actual** strict-core
bridge/export-map object discharging:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_target_v1 (T36/N301),
```

under the discharge acceptance tests of `T148`.

## Background

The lane historically exported:

1. an explicit map-layer nonexport boundary (`T35/F189/P280/N300`),
2. an explicit future-only target naming the missing bridge/export-map object
   (`T36/F190/P281/N301`),
3. stronger post-witness support on the provider-object witness track via an
   actual object-support carrier (`T146/F301/N413`).

On the current repo state, the lane additionally exports an **actual**
strict-core bridge/export-map object (`F311/N422`).

`P388` checks whether the current exports satisfy the discharge acceptance
tests of `T148` and therefore discharge `N301`, while preserving the hard
limits (no theta-source export; `QW-2191` stays open).

## Probe table (T148 acceptance tests)

| Acceptance test (T148) | Verdict | Evidence |
|---|---|---|
| actual map object `Upsilon_*_export_map_object_v1` exported | YES | `F311/N422` export `Upsilon_residual_datum_sigma_int_bridge_export_map_object_v1` as an actual strict-core bridge/export-map object satisfying `T148` |
| typed map-shape `sigma_int_strict_derived_v1 -> residual_orientation_datum_target_slot` exported | YES | `F311` exports `E_sigma_int_to_residual_datum_bridge_export_map_object_v1 : sigma_int_strict_derived_v1 -> residual_orientation_datum_target_slot` (strict source-upgraded sigma-int; no theta inputs) |
| sigma-int strict derivation/source upgrade discharged | YES | `F307/N418` export a strict sigma-int source upgrade on a declared domain (explicit premise-based provenance; no hybrid FR reuse) |
| sigma-int gauge-quotient safety discharged | YES | `F308/N419` export a declared gauge action and an explicit invariance/quotient-level witness (no gauge fixing) |
| selector-track identification beyond overlay-only discharged | YES | `F310/N421` export a strict-core selector-track identification witness (keeps `QW-2191` open) |
| noncyclic contract (no theta/populated inputs) satisfied | YES | `F311/N422` explicitly forbid `theta_{1,2}` and populated basis-pair inputs (respects `N18`) |
| observer-free contract satisfied | YES | `F306..F311` exports are observer-free (no `K_obs` indexing) |
| no implied selector closure / `QW-2191` discipline | YES | `F310/N421` + `F311/N422` explicitly keep `QW-2191` open (no implied `S_sel_int` nor selector closure) |

## Exact verdict

The strict-core bridge/export-map object is now exported:

```text
Upsilon_residual_datum_sigma_int_bridge_export_map_object_target_v1: DISCHARGED (F311/N422)
N300 map-layer nonexport boundary: superseded as a current-state description (historical as-of 2026-03-09)
```

Hard limits remain in force (explicit in `F311/N422`):

1. the exported **map object** remains sign-only: it exports no strict-core theta-source and does not populate the
   residual target slot **by itself**,
2. no implied selector closure; `QW-2191` remains open.

Update (`2026-03-15`, strict-core hygiene; no false pass):

1. Strict-core theta supply on the sigma-int corridor lane is now exported in a slot-free construction class (`T162`)
   satisfying `T159` (`F451/N489`).
2. An audited inhabitant instance populating the `R1` residual orientation datum target slot constructed from that
   slot-free theta-pair source is exported (`P451`).

This does **not** upgrade the strict-core export-map object (`F311/N422`), and it does **not** discharge the missing
post-`T148` actual object-support layer above the exported map object (update: the object-support target `T130/N395` is
now discharged by `F452/N490`; the export-map object itself remains sign-only and is not silently upgraded).

## Consequence (next honest step)

The next honest move is not to relabel the entire residual-datum lane as
closed.

It is to keep the post-`T148` missing layer explicit:

1. theorem-level discharge of the conditional bridge theorem `T2` (beyond probe-level computability), and
2. continued `QW-2191` nonclosure discipline unless a new strict internal selector source is separately exported.
