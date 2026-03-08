# P235 Current Actual Source Topology Selector Witness Probe

Status: `P235_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P235` tests whether the current repo really exports one actual source-side
selector witness, while keeping the result:

1. below basis-independent selector promotion,
2. below quotient-safe `QW-2191` resolution,
3. below current selector closure,
4. below current global `QW-2191` discharge.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Witness under test:

```text
Pi_sel_src_actual_witness_v1 :
tau_src_candidate_v1 -> Sigma_sel_src_target_v1
```

## What P235 checks

`P235` checks only:

1. the witness is explicitly exported,
2. the domain is exactly `tau_src_candidate_v1`,
3. the codomain is exactly `Sigma_sel_src_target_v1`,
4. the selector datum is realized only through the already exported admissible
   preobserver chart lane,
5. observer remains downstream only,
6. no identification `tau_src_candidate_v1 = S_preLM_strict_core_source_object_v1`
   is claimed,
7. the result remains below basis-independent selector promotion,
8. the result remains below quotient-safe `QW-2191` resolution and below
   current selector closure.

## Result

`P235` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_BELOW_BASIS_INDEPENDENCE_AND_QW2191_AFTER_P235
```

This means:

1. the current repo exports one actual source-side selector witness for
   `tau_src_candidate_v1`,
2. the current repo no longer keeps the selector layer only at future-target
   status,
3. but it still does not export basis-independent selector promotion,
   quotient-safe `QW-2191` resolution, or selector closure.

## Hard limits

`P235` does not establish:

1. basis-independent selector promotion,
2. quotient-safe `QW-2191` resolution,
3. current selector closure,
4. current global `QW-2191` discharge,
5. ToE closure.
