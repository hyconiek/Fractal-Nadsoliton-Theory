# P229 Current Actual Source Topology Barrier-Protected Sign Witness Probe

Status: `P229_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P229` tests whether the current repo really exports one actual source-side
barrier-protected sign witness for `tau_src_candidate_v1`, while keeping the
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
Psi_src_barrier_sign_actual_witness_v1 :
tau_src_candidate_v1 -> barrier_protected_sign_class_v1
```

## What P229 checks

`P229` checks only:

1. the packet is explicitly exported,
2. the codomain is exactly `barrier_protected_sign_class_v1`,
3. the already exported positive barrier margin remains present,
4. the already exported local sign-stability witness remains present,
5. the lift stays source-side,
6. observer remains downstream only,
7. the result remains below full source-topology nontriviality,
8. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P229` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P229
```

This means:

1. the current repo exports one actual barrier-protected sign witness for
   `tau_src_candidate_v1`,
2. the current repo no longer keeps the sign layer only at future-subtarget or
   branch-local status,
3. but it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P229` does not establish:

1. full source-topology nontriviality,
2. basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.
