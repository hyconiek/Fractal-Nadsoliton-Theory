# P390 Current Strict T144+T149 Pre-Bridge Discharge Status Probe

Status: `P390_EXECUTED_CURRENT_STRICT_T144_T149_PREBRIDGE_DISCHARGE_STATUS_PROBE_NO_FALSE_PASS`  
As of: `2026-03-10`

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
| `T144` discharged by an actual strict equipartition witness | NO | `P384` records `Omega_16_v1 / mu_eq_v1 / four-bit witness` absent |
| `T149` discharged by an actual strict FR-sign derivation/source-upgrade | NO | `P389` records `C_v1 / pi_1(C_v1) ≅ Z_2 / chi_FR_strict_v1` absent |
| strict bridge/export-map object (`T148`) therefore admissible to claim | NO | `P388` records map object absent; upstream strict prerequisites still missing |

## Exact verdict

```text
T144: NOT DISCHARGED
T149: NOT DISCHARGED
T148 (map export): remains blocked upstream of these discharges
```

## Consequence (next honest step)

The next honest move is not to relabel the state as closed.

It is to discharge at least one strict prerequisite package, e.g.:

1. the strict `C_v1 + pi_1(C_v1) ≅ Z_2 + chi_FR_strict_v1` package (`T149`), or
2. the strict `Omega_16_v1 + mu_eq_v1 + four-bit witness` package (`T144`),

without importing hybrid FR support or legacy-only operator decompositions as
strict sources.

