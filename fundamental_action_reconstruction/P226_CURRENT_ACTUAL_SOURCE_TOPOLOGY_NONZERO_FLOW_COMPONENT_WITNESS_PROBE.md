# P226 Current Actual Source Topology Nonzero-Flow Component Witness Probe

Status: `P226_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P226` tests whether the current repo really exports one actual source-side
scalar nonzero-flow component witness, while keeping the result:

1. below barrier-protected sign discharge,
2. below full source-topology nontriviality discharge,
3. below basis-independent selector promotion,
4. below quotient-safe `QW-2191` resolution,
5. below current selector closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Witness under test:

```text
xi_src_nonzero_flow_component_witness_v1 := |cos(phi)|
```

## What P226 checks

`P226` checks only:

1. the packet is explicitly exported,
2. the witness is positive,
3. the witness stays source-side,
4. observer remains downstream only,
5. the result remains below barrier-sign discharge,
6. the result remains below full source-topology nontriviality,
7. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P226` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_NONZERO_FLOW_COMPONENT_WITNESS_BELOW_BARRIER_SIGN_AND_FULL_NONTRIVIALITY_AFTER_P226
```

This means:

1. the current repo exports one actual scalar nonzero-flow component witness for
   `tau_src_candidate_v1`,
2. but it still does not export barrier-protected sign,
3. and it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P226` does not establish:

1. barrier-protected sign,
2. full source-topology nontriviality,
3. basis-independent selector promotion,
4. quotient-safe `QW-2191` resolution,
5. current selector closure,
6. current global `QW-2191` discharge,
7. ToE closure.
