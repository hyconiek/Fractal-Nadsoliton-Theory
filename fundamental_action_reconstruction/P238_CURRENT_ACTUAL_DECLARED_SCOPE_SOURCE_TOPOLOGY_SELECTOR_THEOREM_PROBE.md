# P238 Current Actual Declared-Scope Source Topology Selector Theorem Probe

Status: `P238_EXECUTED_CURRENT_ACTUAL_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P238` tests whether the current repo really exports one actual declared-scope
Source Topology Selector theorem, while keeping the result:

1. below current strict-core selector closure,
2. below current global selector closure,
3. below current global `QW-2191` discharge,
4. below ToE closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Theorem witness under test:

```text
T14_src_selector_declared_scope_actual_witness_v1 :
tau_src_candidate_v1 -> declared_scope_source_topology_selector_theorem_target_v1
```

## What P238 checks

`P238` checks only:

1. the theorem witness is explicitly exported,
2. the domain is exactly `tau_src_candidate_v1`,
3. actual `T14-L1` full nontriviality support is present,
4. actual `T14-L2` observer-free upstream scope support is present,
5. actual `T14-L3` basis-independent selector-promotion support is present,
6. actual `T14-L4` quotient-safe `QW-2191` support is present,
7. `T14-L5` observer-downstream-only boundary remains explicit,
8. no identification `tau_src_candidate_v1 = S_preLM_strict_core_source_object_v1`
   is claimed,
9. the theorem remains declared-scope only,
10. the result remains below current selector closure and below current global
    `QW-2191` discharge.

## Result

`P238` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_BELOW_CURRENT_SELECTOR_CLOSURE_AFTER_P238
```

This means:

1. the current repo exports one actual declared-scope Source Topology Selector
   theorem witness for `tau_src_candidate_v1`,
2. the `T14` route no longer remains only theorem-spec plus intermediate
   witness layers,
3. but the result still does not export strict-core selector closure,
   global selector closure, or global `QW-2191` discharge.

## Hard limits

`P238` does not establish:

1. current strict-core selector closure,
2. current global selector closure,
3. current global `QW-2191` discharge,
4. ToE closure.
