# P281 Current Residual Datum Sigma-Int Bridge Export Map Object Target Probe

Status: `P281_EXECUTED_CURRENT_RESIDUAL_DATUM_SIGMA_INT_BRIDGE_EXPORT_MAP_OBJECT_TARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Probe whether the residual-datum / `sigma_int_candidate` route now supports
one explicit future-only target object for the missing bridge/export map.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| missing map object is sharply named | YES | `P4/P5` localize the exact missing object |
| source object exists | YES | `B4` exports `sigma_int_candidate` |
| residual codomain slot exists | YES | `C37/R1` keep the target slot in scope |
| candidate-fit exists | YES | `C37` keeps overlay candidate-fit in scope |
| conditional bridge theorem spec exists | YES | `T2` packet is present |
| future object carrier grammar exists | YES | `C40-C46` provide schema/carrier grammar |
| actual bridge/export map exists | NO | still absent |
| future-only target object exists | YES | strongest honest positive result below actual export |

## Probe result

`P281` returns:

```text
future-only bridge/export-map object target present: yes
actual bridge/export map present: no
```

## Consequence

The strongest honest current repo reading is:

1. the route now carries one explicit future-only target object for the
   missing map layer,
2. but the exact map itself is still unexported,
3. so the route remains below actual bridge/export-map discharge.
