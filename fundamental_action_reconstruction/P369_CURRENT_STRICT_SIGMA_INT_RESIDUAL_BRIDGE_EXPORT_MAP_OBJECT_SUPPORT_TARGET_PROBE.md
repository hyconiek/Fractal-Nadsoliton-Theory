# P369 Current Strict Sigma-Int Residual Bridge/Export-Map Actual Object Support Target Probe

Status: `P369_EXECUTED_CURRENT_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_TARGET_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Probe whether the current repo exports:

```text
actual bridge/export-map object support (above the exported map object)
```

or only the weaker, admissible result:

```text
one explicit future-only target object naming the missing post-witness
actual object-support layer (F283),
with explicit acceptance tests (T130).
```

On the updated repo state, the strict sigma-int lane already exports an
**actual** strict-core bridge/export-map object satisfying `T148` (`F311/N422`),
so the only remaining question here is the missing **object-support layer**
above that map object (`N302`).

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| export-map object present | YES | `F311/N422` |
| object-support projection layer present | YES | `N385` |
| object-support witness present | YES | `N387` |
| actual object support present | NO | `N302` remains in force (post-`T148` reading: support above the exported map object is absent) |
| actual object-support target named | YES | `T130/F283` |

## Exact verdict

The strongest honest current verdict is:

```text
actual bridge/export-map object support: absent
future-only object-support target naming: exported
```
