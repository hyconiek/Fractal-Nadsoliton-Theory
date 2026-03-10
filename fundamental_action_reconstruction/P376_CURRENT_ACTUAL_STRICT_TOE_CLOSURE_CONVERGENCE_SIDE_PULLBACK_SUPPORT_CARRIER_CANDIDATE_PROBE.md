# P376 Current Actual Strict ToE Closure Convergence-Side Pullback Support Carrier Candidate Probe

Status: `P376_EXECUTED_CURRENT_ACTUAL_STRICT_TOE_CLOSURE_CONVERGENCE_SIDE_PULLBACK_SUPPORT_CARRIER_CANDIDATE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Probe whether the current repo exports one explicit convergence-side pullback
support carrier **candidate** packaging:

1. the provider-object carrier projection lane (`N399`),
2. the sigma-int low-delta projection lane (`N400`),

while remaining strictly below:

- discharge of `N302`,
- any export-map object export,
- any selector closure / `QW-2191` discharge.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| provider positive-window projection candidate exported | YES | `N399` |
| sigma-int positive-window low-delta projection candidate exported | YES | `N400` |
| joint coherence witness candidate exported | YES | `N401` |
| pullback support carrier candidate exported | YES | `T135/F290` |
| actual bridge/export-map object support exported | NO | `N302` remains in force |

## Exact verdict

The strongest honest current verdict is:

```text
pullback support carrier candidate present: YES
actual residual bridge/export-map object support: NO
```

