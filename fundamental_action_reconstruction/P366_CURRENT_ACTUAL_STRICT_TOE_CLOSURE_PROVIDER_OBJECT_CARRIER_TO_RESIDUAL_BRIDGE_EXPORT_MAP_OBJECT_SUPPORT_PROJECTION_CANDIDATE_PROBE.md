# P366 Current Actual Strict ToE Closure Provider-Object Carrier to Residual Bridge/Export-Map Object-Support Projection Candidate Probe

Status: `P366_EXECUTED_CURRENT_ACTUAL_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_TO_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Probe whether the current repo already exports:

```text
an actual provider-carrier -> residual-bridge object-support lift
```

or only the weaker result:

```text
an explicit bridge-facing projection *candidate* (F280)
that remains below actual object support (N302).
```

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| provider-object carrier candidate exported | YES | `N391` |
| strict phasor reduction candidate exported | YES | `N380` |
| bridge-facing projection candidate exported | YES | `T127/F280` |
| actual bridge/export-map object support exported | NO | `N302` remains in force |
| export-map object exported | YES | `F311/N422` discharge `T148` (`N300` superseded; `N301` discharged) |

## Exact verdict

The strongest honest current verdict is:

```text
projection candidate present: YES
actual residual bridge/export-map object support: NO
```
