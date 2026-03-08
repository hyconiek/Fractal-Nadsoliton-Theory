# P217 Current Source Topology Invariant Nontriviality Target Probe

Status: `P217_EXECUTED_CURRENT_SOURCE_TOPOLOGY_INVARIANT_NONTRIVIALITY_TARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P217` tests whether the current repo really exports the future-only
nontriviality target packet introduced by `F129`, while keeping the result:

1. below actual invariant discharge,
2. below selector promotion,
3. below quotient-safe `QW-2191` resolution,
4. below current selector closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Target under test:

```text
Nu_src_nontriv_target_v1 : tau_src_candidate_v1 -> Lambda_src_nontriv_target_v1
```

with:

```text
Lambda_src_nontriv_target_v1 :=
(
  source_limit_nonzero_flow_class_v1,
  barrier_protected_sign_class_v1,
  observer_free_scope_tag_v1
)
```

## What P217 checks

`P217` checks only:

1. the packet is explicitly exported,
2. the domain is `tau_src_candidate_v1`,
3. the codomain is an observer-free source-side nontriviality witness class,
4. the result remains future-only,
5. the result remains below actual nontriviality discharge,
6. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P217` returns:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_INVARIANT_NONTRIVIALITY_TARGET_BELOW_ACTUAL_NONTRIVIALITY_DISCHARGE_AFTER_P217
```

This means:

1. the current repo exports one explicit future-only nontriviality target for
   `tau_src_candidate_v1`,
2. but it still does not export a discharged non-trivial source-topology
   invariant,
3. and it still does not export selector promotion or `QW-2191` discharge.

## Hard limits

`P217` does not establish:

1. actual non-triviality of `tau_src_candidate_v1`,
2. basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.
