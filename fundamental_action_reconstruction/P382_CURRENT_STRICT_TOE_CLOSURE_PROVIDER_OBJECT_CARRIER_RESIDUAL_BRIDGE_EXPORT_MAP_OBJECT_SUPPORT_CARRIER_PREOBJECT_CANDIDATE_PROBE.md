# P382 Current Strict ToE Closure Provider-Object Carrier Residual Bridge/Export-Map Object-Support Carrier Preobject Candidate Probe

Status: `P382_EXECUTED_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_CARRIER_PREOBJECT_CANDIDATE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Probe whether the current repo exports:

1. any **actual** bridge/export-map object-support carrier above the
   provider-object residual witness layer (`N405`),

or only the weaker admissible result:

2. one explicit **candidate preobject** construction for such a carrier
   (`F296/T141`),
   while keeping `N302` in force.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| provider-object carrier candidate present | YES | `N391` |
| provider-object → residual projection layer present | YES | `N398` |
| provider-object residual object-support witness present | YES | `N405` |
| post-witness carrier target named | YES | `N407` |
| actual post-witness object-support carrier present | NO | `N302` remains in force |
| preobject carrier candidate exported | YES | `T141/F296` |

## Exact verdict

The strongest honest current verdict is:

```text
actual post-witness object-support carrier: absent
explicit preobject carrier candidate: exported
```

