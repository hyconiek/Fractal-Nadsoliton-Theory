# P374 Current Actual Strict Sigma-Int To Residual Datum Bridge Export-Map Object-Support Projection Candidate (Positive-Window Low-Delta) Probe

Status: `P374_EXECUTED_CURRENT_ACTUAL_STRICT_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_POSITIVE_WINDOW_LOW_DELTA_PROBE_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Probe whether the current repo exports a **low-delta** corridor-protected
sigma-int projection candidate suitable for convergence-side coherence
comparison, while remaining explicitly below:

1. discharge of `N302`,
2. bridge/export-map export,
3. strict-core theta export,
4. selector closure / `QW-2191` discharge.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| base positive-window sigma-int projection candidate exported | YES | `N384` |
| low-delta positive-window sigma-int projection candidate exported | YES | `T133/F288` |
| actual bridge/export-map object support exported | NO | `N302` remains in force |
| export-map object exported | NO | `N300/N301` scope preserved |

## Exact verdict

The strongest honest current verdict is:

```text
low-delta positive-window sigma-int projection candidate present: YES
actual residual bridge/export-map object support: NO
```

