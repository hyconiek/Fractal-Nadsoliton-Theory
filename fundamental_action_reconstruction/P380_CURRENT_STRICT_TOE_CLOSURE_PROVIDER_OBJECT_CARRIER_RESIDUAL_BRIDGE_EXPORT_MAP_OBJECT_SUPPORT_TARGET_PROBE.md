# P380 Current Strict ToE Closure Provider-Object Carrier Residual Bridge/Export-Map Actual Object Support Target Probe

Status: `P380_EXECUTED_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_TARGET_PROBE_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Probe whether the current provider-object carrier residual route exports:

```text
actual bridge/export-map object support
```

or only the weaker, admissible result:

```text
one explicit future-only target object naming the missing post-witness
actual object-support layer (F294),
with explicit acceptance tests (T139).
```

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| object-support projection layer present | YES | `N398` |
| object-support witness present | YES | `N405` |
| actual post-witness object-support carrier present | YES | `T146/F301/N413` exports a carrier above `N405` (below `N300`) |
| actual object-support target named | YES | `T139/F294` |

## Exact verdict

The strongest honest current verdict is:

```text
post-witness object-support carrier: present (provider-object witness track)
object-support target naming: exported
```
