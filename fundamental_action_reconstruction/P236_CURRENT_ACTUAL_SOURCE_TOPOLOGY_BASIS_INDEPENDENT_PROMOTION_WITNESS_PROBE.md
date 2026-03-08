# P236 Current Actual Source Topology Basis Independent Promotion Witness Probe

Status: `P236_EXECUTED_CURRENT_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P236` tests whether the current repo really exports the actual
basis-independent selector-promotion witness introduced by `F148`, while
keeping the result:

1. below quotient-safe `QW-2191` resolution,
2. below current selector closure,
3. below current global `QW-2191` discharge.

## Fixed input

Input witness:

```text
Upsilon_sel_basis_actual_witness_v1 :
tau_src_candidate_v1 -> Sigma_sel_basis_free_target_v1
```

## What P236 checks

`P236` checks only:

1. the actual basis-independent promotion witness is explicitly exported,
2. the input packet remains `tau_src_candidate_v1`,
3. the codomain packet remains `Sigma_sel_basis_free_target_v1`,
4. the basis-free axis class is exported,
5. the basis-free signed-split class is exported,
6. the basis-free preobserver scope tag is exported,
7. the witness remains source-side and observer-free in the witness domain,
8. the witness remains kernel-split-safe,
9. the witness does not identify `tau_src_candidate_v1` with
   `S_preLM_strict_core_source_object_v1`,
10. the result remains below quotient-safe `QW-2191` resolution and below
    current selector closure.

## Result

`P236` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_BELOW_QW2191_AFTER_P236
```

This means:

1. the current repo exports one actual basis-independent selector-promotion
   witness from `tau_src_candidate_v1` into a basis-free selector packet,
2. but it still does not export quotient-safe `QW-2191` resolution,
3. and it still does not export current selector closure.

## Hard limits

`P236` does not establish:

1. quotient-safe `QW-2191` resolution,
2. current selector closure,
3. current global `QW-2191` discharge,
4. ToE closure.
