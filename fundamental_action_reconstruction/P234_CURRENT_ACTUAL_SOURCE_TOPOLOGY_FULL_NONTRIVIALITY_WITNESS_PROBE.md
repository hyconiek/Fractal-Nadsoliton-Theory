# P234 Current Actual Source Topology Full Nontriviality Witness Probe

Status: `P234_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P234` tests whether the current repo really exports one actual full
source-topology nontriviality witness, while keeping the result:

1. below actual source-side selector witness,
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
Theta_src_nontriv_actual_discharge_witness_v1 :
tau_src_candidate_v1 -> actual_full_source_topology_nontriviality_discharge_target_v1
```

## What P234 checks

`P234` checks only:

1. the witness is explicitly exported,
2. the domain is exactly `tau_src_candidate_v1`,
3. the codomain is exactly
   `actual_full_source_topology_nontriviality_discharge_target_v1`,
4. the already exported actual assembly witness remains present,
5. the witness stays source-side,
6. observer remains downstream only,
7. the result remains below actual source-side selector witness,
8. the result remains below basis-independent selector promotion and below
   `QW-2191` discharge.

## Result

`P234` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_WITNESS_BELOW_SELECTOR_WITNESS_AND_BASIS_INDEPENDENCE_AFTER_P234
```

This means:

1. the current repo exports one actual full source-topology nontriviality
   witness for `tau_src_candidate_v1`,
2. the current repo no longer keeps the full nontriviality layer only at
   future-target or assembly-only status,
3. but it still does not export an actual source-side selector witness,
   basis-independent selector promotion, or `QW-2191` discharge.

## Hard limits

`P234` does not establish:

1. actual source-side selector witness,
2. basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.
