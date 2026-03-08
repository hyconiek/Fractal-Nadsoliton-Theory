# P231 Current Actual Source Topology Nonzero-Flow Witness Probe

Status: `P231_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P231` tests whether the current repo really exports one actual source-side
nonzero-flow witness for `tau_src_candidate_v1`, while keeping the result:

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
Xi_src_nonzero_flow_actual_witness_v1 :
tau_src_candidate_v1 -> source_limit_nonzero_flow_class_v1
```

## What P231 checks

`P231` checks only:

1. the packet is explicitly exported,
2. the codomain is exactly `source_limit_nonzero_flow_class_v1`,
3. the already exported scalar component witness remains positive,
4. the already exported source-limit support packet remains present,
5. the lift stays source-side,
6. observer remains downstream only,
7. the result remains below full source-topology nontriviality,
8. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P231` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P231
```

This means:

1. the current repo exports one actual nonzero-flow witness for
   `tau_src_candidate_v1`,
2. the current repo no longer keeps the nonzero-flow layer only at
   future-subtarget or scalar-component status,
3. but it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P231` does not establish:

1. full source-topology nontriviality,
2. basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.
