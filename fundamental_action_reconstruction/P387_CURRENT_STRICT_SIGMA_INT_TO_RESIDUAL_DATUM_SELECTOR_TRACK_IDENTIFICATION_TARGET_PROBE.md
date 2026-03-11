# P387 Current Strict Sigma-Int to Residual Datum Selector-Track Identification Target Probe

Status: `P387_EXECUTED_CURRENT_STRICT_SIGMA_INT_TO_RESIDUAL_DATUM_SELECTOR_TRACK_IDENTIFICATION_TARGET_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `R1/P5/N8`, the strict-core residual-datum bridge lane remained blocked
by a finite list of missing objects.

On the current repo state, the following previously-missing prerequisite
packages are now exported:

1. sigma-int strict derivation/source upgrade (`F307/N418`),
2. sigma-int gauge-quotient safety (`F308/N419`),
3. bridge/export-map object actual map object export (`F311/N422`).

`P5` still lists one additional missing ingredient:

```text
selector-track identification beyond overlay-only compatibility
```

`P387` probes whether the current repo already exports any strict-core witness
of that kind, or whether the strongest honest positive move is only target
naming (`T147`).

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| residual `Z2` slot separated and named | YES | `B6/C37` |
| candidate-fit `sigma_int_candidate ~ residual Z2 slot` present | YES | `C37` |
| selector-track compatibility present only as overlay | YES | `B7` + `N8` |
| strict-core selector-track identification witness exported | YES | `F310/N421` export `Chi_sigma_int_residual_datum_selector_track_identification_witness_v1` |
| future-only target naming still the strongest positive move | NO | the missing ingredient is no longer “only nameable”; an actual strict witness is now exported (`N421`) |

## Exact verdict

The strongest honest current verdict is:

```text
strict-core selector-track identification beyond overlay-only compatibility: present (QW-2191 remains open)
future-only target naming: superseded by actual export (still admissible as a target name)
```

## Finite missing-object list (selector-track identification sublane)

After `F310/N421`, this sublane no longer lacks a strict-core selector-track
identification witness of the form required by `T147`.

The remaining missing objects live **outside** this sublane, e.g.:

1. strict theta-source export and/or residual target-slot population (still
   absent),
2. no implied selector closure / `QW-2191` discharge discipline (must remain
   explicit unless separately discharged).

## Consequence (next honest step)

If this sublane is continued, the next honest move is:

1. keep the residual bridge lane from being misread as theta/population export
   (those layers remain absent),
2. keep `QW-2191` explicitly open unless a new strict selector/symmetry-breaking
   ingredient is separately exported.
