# P220 Current Source Topology Observer-Free Scope Subtarget Probe

Status: `P220_EXECUTED_CURRENT_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P220` tests whether the current repo really exports the future-only
observer-free scope subtarget packet introduced by `F132`, while keeping the
result:

1. below actual observer-free scope discharge,
2. below actual barrier-protected sign discharge,
3. below actual nonzero-flow discharge,
4. below full source-topology nontriviality,
5. below selector promotion,
6. below quotient-safe `QW-2191` resolution,
7. below current selector closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Subtarget under test:

```text
Omega_src_observer_free_scope_target_v1 : tau_src_candidate_v1 -> observer_free_scope_tag_v1
```

## What P220 checks

`P220` checks only:

1. the packet is explicitly exported,
2. the domain is `tau_src_candidate_v1`,
3. the codomain is `observer_free_scope_tag_v1`,
4. the result remains observer-free,
5. the result remains future-only,
6. the result remains below actual observer-free scope discharge,
7. the result remains below actual barrier-protected sign discharge,
8. the result remains below actual nonzero-flow discharge,
9. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P220` returns:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_SUBTARGET_BELOW_ACTUAL_OBSERVER_FREE_SCOPE_DISCHARGE_AFTER_P220
```

This means:

1. the current repo exports one explicit future-only observer-free scope
   subtarget for `tau_src_candidate_v1`,
2. but it still does not export an actual discharged observer-free scope
   witness,
3. and it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P220` does not establish:

1. actual observer-free scope of `tau_src_candidate_v1`,
2. actual barrier-protected sign of `tau_src_candidate_v1`,
3. actual nonzero-flow of `tau_src_candidate_v1`,
4. full source-topology nontriviality,
5. basis-independent selector promotion,
6. quotient-safe `QW-2191` resolution,
7. current selector closure,
8. current global `QW-2191` discharge,
9. ToE closure.
