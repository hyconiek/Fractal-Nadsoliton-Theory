# P264 Current Strict Source-To-Pair-Population Noncyclic Anchor Target Probe

Status: `P264_EXECUTED_CURRENT_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P264` tests whether the current repo really exports one explicit future-only
strict-side solution target of the form described in `T26/F173`, while
keeping the result:

1. below actual strict-core ingredient discharge,
2. below actual `E_orient`,
3. below admissible `S_sel_int`,
4. below strict-core selector closure,
5. below ToE closure.

## Fixed target under test

```text
V_strict_src_pair_population_anchor_target_v1 :
  tau_src_candidate_v1 -> strict_src_pair_population_noncyclic_anchor_target_v1
```

## What P264 checks

`P264` checks only:

1. the target is explicitly exported as future-only,
2. the target is source-side rather than observer-side,
3. the target does not use `K_obs` as a primary selector source,
4. the target does not reactivate `K_legacy_ont` as the live forward
   constructive kernel,
5. the target is pair-indexed rather than only class-level,
6. the target is explicitly motivated as a proposed noncyclic anchor against
   the `N18` loop,
7. the target is explicitly motivated as a proposed new ingredient class
   against the `N283` official-lane freeze,
8. the target remains below actual theta supply, actual pair population,
   actual `E_orient`, admissible `S_sel_int`, and closure.

## Result

`P264` returns:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_AFTER_P264
```

This means:

1. the current repo now contains one explicit rigorous proposal for the
   missing strict-side ingredient class,
2. the proposal is sharper than a generic "new ingredient needed" statement,
3. the proposal still remains only a target,
4. no actual strict-core ingredient is yet exported.

## Hard limits

`P264` does not establish:

1. actual source-side anchor discharge,
2. actual theta values,
3. actual populated basis-pair instance,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. actual ToE closure.
