# P230 Current Actual Source Topology Observer-Free Scope Witness Probe

Status: `P230_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P230` tests whether the current repo really exports one actual source-side
observer-free scope witness for `tau_src_candidate_v1`, while keeping the
result:

1. below full source-topology nontriviality discharge,
2. below basis-independent selector promotion,
3. below quotient-safe `QW-2191` resolution,
4. below current selector closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Witness under test:

```text
Omega_src_observer_free_scope_actual_witness_v1 :
tau_src_candidate_v1 -> observer_free_scope_tag_v1
```

## What P230 checks

`P230` checks only:

1. the packet is explicitly exported,
2. the codomain is exactly `observer_free_scope_tag_v1`,
3. the source packet remains upstream of observer,
4. the already exported observer-to-upstream nonparticipation support remains
   present,
5. the theorem-level downstream-only observer support from `N163/N234`
   remains present,
6. the lift stays source-side,
7. the result remains below full source-topology nontriviality,
8. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P230` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P230
```

This means:

1. the current repo exports one actual observer-free scope witness for
   `tau_src_candidate_v1`,
2. the current repo no longer keeps the scope layer only at future-subtarget
   status,
3. but it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P230` does not establish:

1. full source-topology nontriviality,
2. basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.
