# P219 Current Source Topology Barrier-Protected Sign Subtarget Probe

Status: `P219_EXECUTED_CURRENT_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_SUBTARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P219` tests whether the current repo really exports the future-only
barrier-protected sign subtarget packet introduced by `F131`, while keeping the
result:

1. below actual barrier-protected sign discharge,
2. below actual nonzero-flow discharge,
3. below full source-topology nontriviality,
4. below selector promotion,
5. below quotient-safe `QW-2191` resolution,
6. below current selector closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Subtarget under test:

```text
Psi_src_barrier_sign_target_v1 : tau_src_candidate_v1 -> barrier_protected_sign_class_v1
```

## What P219 checks

`P219` checks only:

1. the packet is explicitly exported,
2. the domain is `tau_src_candidate_v1`,
3. the codomain is `barrier_protected_sign_class_v1`,
4. the result remains observer-free,
5. the result remains future-only,
6. the result remains below actual barrier-protected sign discharge,
7. the result remains below actual nonzero-flow discharge,
8. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P219` returns:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_SUBTARGET_BELOW_ACTUAL_BARRIER_PROTECTED_SIGN_DISCHARGE_AFTER_P219
```

This means:

1. the current repo exports one explicit future-only barrier-protected sign
   subtarget for `tau_src_candidate_v1`,
2. but it still does not export an actual discharged barrier-protected sign
   witness,
3. and it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P219` does not establish:

1. actual barrier-protected sign of `tau_src_candidate_v1`,
2. actual nonzero-flow of `tau_src_candidate_v1`,
3. full source-topology nontriviality,
4. basis-independent selector promotion,
5. quotient-safe `QW-2191` resolution,
6. current selector closure,
7. current global `QW-2191` discharge,
8. ToE closure.
