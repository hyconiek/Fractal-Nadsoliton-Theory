# P390 Current Strict T144+T149 Pre-Bridge Discharge Status Probe

Status: `P390_EXECUTED_CURRENT_STRICT_T144_T149_PREBRIDGE_DISCHARGE_STATUS_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Provide one combined probe confirming whether the current repo state already
exports **actual** strict discharges of:

1. `T144` (strict derivation/source-upgrade of `alpha_geo = 4 ln 2` via a
   16-microstate equipartition witness), and
2. `T149` (strict derivation/source-upgrade of the sigma-int sign candidate
   `chi_FR(gamma_pi1) ∈ {+1,-1}` on a declared strict domain).

This probe exists to prevent “roadmap wording” from being mistaken for
discharge.

## Probe table

| Check | Verdict | Evidence |
|---|---|---|
| `T144` discharged by an actual strict equipartition witness | YES | `F309/N420` export `Omega_16_v1`, a symmetry-forced equipartition measure, a four-bit witness, and `alpha_geo_strict_derived_v1 := H(mu_eq_v1) = 4 ln 2` |
| `T149` discharged by an actual strict FR-sign source-upgrade package | YES | `F307/N418` export `chi_FR_strict_v1` (explicit strict-side premise, no hybrid reuse) and `sigma_int_strict_derived_v1` on a strict domain with `pi_1(C_v1) ≅ Z_2` |
| strict bridge/export-map object (`T148`) discharged as an actual export-map object | YES | `F311/N422` export `Upsilon_residual_datum_sigma_int_bridge_export_map_object_v1` satisfying `T148` (residual `Z2` population only; no theta export; `QW-2191` remains open) |

## Exact verdict

```text
T144: DISCHARGED (symmetry-forced equipartition witness; strict alpha-geo source upgrade exported)
T149: DISCHARGED (premise-based strict source upgrade; no hybrid FR reuse)
T148 (map export): DISCHARGED (actual strict-core export-map object exported; residual Z2 population only)
```

## Consequence (next honest step)

The next honest move is not to relabel the state as closed.

It is to keep the post-`T148` missing layers explicit, e.g.:

1. strict theta-source export and/or residual target-slot population (still
   absent),
2. continued `QW-2191` nonclosure discipline unless a new strict
   selector/symmetry-breaking ingredient is separately exported,

without importing hybrid FR support or legacy-only operator decompositions as
strict sources.
