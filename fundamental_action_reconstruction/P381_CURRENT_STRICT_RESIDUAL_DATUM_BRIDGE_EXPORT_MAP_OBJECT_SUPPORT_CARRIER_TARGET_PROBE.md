# P381 Current Strict Residual Datum Bridge/Export-Map Actual Object-Support Carrier Target Probe

Status: `P381_EXECUTED_CURRENT_STRICT_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_CARRIER_TARGET_PROBE_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Probe whether the current repo exports:

1. any **actual** bridge/export-map object-support carrier above the current
   witness layers (`N387` / `N405`),

or only the weaker admissible result:

2. one explicit future-only target object naming that missing carrier layer
   (`F295` / `T140`),
   while keeping `N302` in force.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| sigma-int object-support witness present | YES | `N387` |
| provider-object object-support witness present | YES | `N405` |
| actual object-support carrier above witness present | NO | `N302` remains in force |
| object-support carrier target named | YES | `T140/F295` |

## Exact verdict

The strongest honest current verdict is:

```text
actual bridge/export-map object-support carrier above witness: absent
future-only object-support carrier target naming: exported
```

