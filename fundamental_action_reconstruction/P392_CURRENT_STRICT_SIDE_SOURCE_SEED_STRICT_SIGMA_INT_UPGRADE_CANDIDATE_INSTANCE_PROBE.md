# P392 Current Strict-Side Source-Seed Strict Sigma-Int Upgrade Candidate Instance Probe

Status: `P392_EXECUTED_CURRENT_STRICT_SIDE_SOURCE_SEED_STRICT_SIGMA_INT_UPGRADE_CANDIDATE_INSTANCE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Probe whether the current repo exports the strict-side seed candidate instance:

```text
S_sel_int_candidate_seed_v1 :=
(QW-2206_local_topological_protection_layer_in_local_B_tilde_1_sector, sigma_int_strict_derived_v1)
```

as a strict-sigma-int-upgraded replacement for the sigma-int slot of the
historical seed `S_sel_int_candidate_seed_v0` (`F36`),
while keeping the result explicitly below admissible `S_sel_int`,
selector closure, and `QW-2191` discharge.

## Probe table

| Check | Verdict | Evidence |
|---|---|---|
| strict sigma-int source-upgrade datum exported | YES | `F307/N418` export `sigma_int_strict_derived_v1` (explicit premise; no hybrid reuse) |
| strict-side seed candidate instance `S_sel_int_candidate_seed_v1` exported | YES | `F318` |
| no implied admissible `S_sel_int` / selector closure | YES | `F318` hard limits keep `S_sel_int` and `QW-2191` explicitly open |

## Result

`P392` returns:

```text
CURRENT_REPO_EXPORTS_ONE_STRICT_SIGMA_INT_UPGRADED_STRICT_SIDE_SOURCE_SEED_CANDIDATE_INSTANCE_AFTER_P392
```

## Hard limits

`P392` does not establish:

1. admissible `S_sel_int`,
2. strict-core selector closure,
3. `QW-2191` discharge,
4. strict-core theta export,
5. ToE closure.

