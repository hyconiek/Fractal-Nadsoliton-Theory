# P227 Current Actual Source Topology Barrier Sign Component Witness Probe

Status: `P227_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P227` tests whether the current repo really exports one actual source-side
scalar barrier-sign component witness, while keeping the result:

1. below full barrier-protected sign discharge,
2. below full source-topology nontriviality discharge,
3. below basis-independent selector promotion,
4. below quotient-safe `QW-2191` resolution,
5. below current selector closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Witnesses under test:

```text
delta_src_barrier_sign_margin_v1 := (pi/2) - |phi|
psi_src_barrier_sign_component_witness_v1 := sign(cos(phi))
```

## What P227 checks

`P227` checks only:

1. the packet is explicitly exported,
2. the already exported nonzero-flow component support remains positive,
3. the barrier margin is positive,
4. the sign witness equals `+1` on the current declared core branch,
5. the witness stays source-side,
6. observer remains downstream only,
7. the result remains below full barrier-protected sign discharge,
8. the result remains below full source-topology nontriviality,
9. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P227` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BARRIER_SIGN_COMPONENT_WITNESS_BELOW_FULL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P227
```

This means:

1. the current repo exports one actual scalar barrier-sign component witness for
   `tau_src_candidate_v1`,
2. the current repo exports one explicit positive margin to the first
   sign-flip boundary on the declared core branch,
3. but it still does not export full barrier-protected sign discharge,
4. and it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P227` does not establish:

1. full barrier-protected sign of `tau_src_candidate_v1`,
2. full source-topology nontriviality,
3. basis-independent selector promotion,
4. quotient-safe `QW-2191` resolution,
5. current selector closure,
6. current global `QW-2191` discharge,
7. ToE closure.
