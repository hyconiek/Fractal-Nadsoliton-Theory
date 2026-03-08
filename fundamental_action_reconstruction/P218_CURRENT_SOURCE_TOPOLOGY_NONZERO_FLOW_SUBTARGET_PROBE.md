# P218 Current Source Topology Nonzero-Flow Subtarget Probe

Status: `P218_EXECUTED_CURRENT_SOURCE_TOPOLOGY_NONZERO_FLOW_SUBTARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P218` tests whether the current repo really exports the future-only nonzero-flow
subtarget packet introduced by `F130`, while keeping the result:

1. below actual nonzero-flow discharge,
2. below full source-topology nontriviality,
3. below selector promotion,
4. below quotient-safe `QW-2191` resolution,
5. below current selector closure.

## Fixed input

Input packet:

```text
tau_src_candidate_v1
```

Subtarget under test:

```text
Xi_src_nonzero_flow_target_v1 : tau_src_candidate_v1 -> source_limit_nonzero_flow_class_v1
```

## What P218 checks

`P218` checks only:

1. the packet is explicitly exported,
2. the domain is `tau_src_candidate_v1`,
3. the codomain is `source_limit_nonzero_flow_class_v1`,
4. the result remains observer-free,
5. the result remains future-only,
6. the result remains below actual nonzero-flow discharge,
7. the result remains below selector promotion and below `QW-2191` discharge.

## Result

`P218` returns:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_NONZERO_FLOW_SUBTARGET_BELOW_ACTUAL_NONZERO_FLOW_DISCHARGE_AFTER_P218
```

This means:

1. the current repo exports one explicit future-only nonzero-flow subtarget for
   `tau_src_candidate_v1`,
2. but it still does not export an actual discharged nonzero-flow source
   invariant,
3. and it still does not export full source-topology nontriviality, selector
   promotion, or `QW-2191` discharge.

## Hard limits

`P218` does not establish:

1. actual nonzero-flow of `tau_src_candidate_v1`,
2. barrier-protected sign,
3. full source-topology nontriviality,
4. basis-independent selector promotion,
5. quotient-safe `QW-2191` resolution,
6. current selector closure,
7. current global `QW-2191` discharge,
8. ToE closure.
