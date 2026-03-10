# P280 Current Actual Residual Datum Sigma-Int Bridge Export Map Nonexport Boundary Probe

Status: `P280_EXECUTED_CURRENT_ACTUAL_RESIDUAL_DATUM_SIGMA_INT_BRIDGE_EXPORT_MAP_NONEXPORT_BOUNDARY_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Probe whether the residual-datum / `sigma_int_candidate` route already exports
the actual bridge/export map itself, or whether the strongest honest reading
remains one exact nonexport boundary on that layer.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| `sigma_int_candidate` exists | YES | `B4` exports the candidate datum |
| residual orientation datum slot exists | YES | `C37/R1` keep the codomain slot in scope |
| candidate-fit exists | YES | `C37` keeps overlay candidate-fit in scope |
| conditional bridge theorem spec exists | YES | `T2` packet is present |
| actual bridge-map layer support exists | YES | `F188/P279/N299` already export target support |
| actual bridge/export map exists | NO | still absent |
| actual theta source exists | NO | still absent |
| exact nonexport boundary on the map layer exists | YES | strongest honest current reading below discharge |

## Probe result

`P280` returns:

```text
actual bridge/export map present: no
exact bridge/export-map nonexport boundary present: yes
```

## Consequence

The strongest honest current repo reading is:

1. the route is already above target-only and above support-free status,
2. but the exact bridge/export map itself is still unexported,
3. so the bridge/export-map layer is now frozen as one route-specific
   current-state nonexport boundary.
