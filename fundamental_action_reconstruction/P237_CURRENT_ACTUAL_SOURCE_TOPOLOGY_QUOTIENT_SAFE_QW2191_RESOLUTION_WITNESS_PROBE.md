# P237 Current Actual Source Topology Quotient-Safe QW2191 Resolution Witness Probe

Status: `P237_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P237` tests whether the current repo really exports the actual
quotient-safe `QW-2191` resolution witness introduced by `F149`, while
keeping the result:

1. below current selector closure,
2. below current global `QW-2191` discharge,
3. below ToE closure.

## Fixed input

Input witness:

```text
Phi_qw2191_safe_actual_witness_v1 :
tau_src_candidate_v1 -> actual_quotient_safe_qw2191_resolution_target_v1
```

## What P237 checks

`P237` checks only:

1. the actual quotient-safe resolution witness is explicitly exported,
2. the input packet remains `tau_src_candidate_v1`,
3. the codomain target remains
   `actual_quotient_safe_qw2191_resolution_target_v1`,
4. the actual basis-independent selector-promotion witness is present,
5. the `QW-2191` kernel-alone obstruction remains explicit,
6. the witness resolves the ambiguity only to a distinguished quotient class,
   not to a raw unique `theta`,
7. the witness remains source-side and observer-free in the witness domain,
8. the witness remains kernel-split-safe,
9. the result remains below current selector closure and below current global
   `QW-2191` discharge.

## Result

`P237` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_WITNESS_BELOW_CURRENT_SELECTOR_CLOSURE_AFTER_P237
```

This means:

1. the current repo exports one actual source-side quotient-safe
   `QW-2191` resolution witness in the declared source-topology scope,
2. but it still does not export current selector closure,
3. and it still does not export current global `QW-2191` discharge.

## Hard limits

`P237` does not establish:

1. current selector closure,
2. current global `QW-2191` discharge,
3. ToE closure.
