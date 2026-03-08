# P224 Current Source Topology Basis-Independent Promotion Target Probe

Status: `P224_EXECUTED_CURRENT_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_TARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P224` tests whether the current repo really exports the future-only
basis-independent promotion target introduced by `F136`, while keeping the
result:

1. below actual basis-independent selector promotion,
2. below quotient-safe `QW-2191` resolution,
3. below current selector closure.

## Fixed input

Input promotion target:

```text
Upsilon_sel_basis_target_v1 :
(Theta_src_nontriv_discharge_target_v1, Pi_sel_src_target_v1)
  -> Sigma_sel_basis_free_target_v1
```

## What P224 checks

`P224` checks only:

1. the basis-independent promotion target is explicitly exported,
2. the domain target `Theta_src_nontriv_discharge_target_v1` is present,
3. the domain target `Pi_sel_src_target_v1` is present,
4. the codomain target `Sigma_sel_basis_free_target_v1` is present,
5. the target remains source-side,
6. the target remains future-only,
7. the target remains below actual basis-independent selector promotion and
   below quotient-safe `QW-2191`.

## Result

`P224` returns:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_TARGET_BELOW_QW2191_QUOTIENT_SAFETY_AFTER_P224
```

This means:

1. the current repo exports one explicit future-only target from the pair
   `(Theta_src_nontriv_discharge_target_v1, Pi_sel_src_target_v1)` to a later
   basis-free selector packet,
2. but it still does not export actual basis-independent selector promotion,
3. and it still does not export quotient-safe `QW-2191` discharge.

## Hard limits

`P224` does not establish:

1. actual full source-topology nontriviality,
2. actual basis-independent selector promotion,
3. quotient-safe `QW-2191` resolution,
4. current selector closure,
5. current global `QW-2191` discharge,
6. ToE closure.
