# P71 Strict-Side Replaced-Branch Probe For Legacy Weinberg Role

Status: `P71_EXECUTED_STRICT_SIDE_REPLACED_BRANCH_PROBE_FOR_LEGACY_WEINBERG_ROLE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N73`, the next honest direct question is:

```text
does the current repo export an explicit strict-side replaced verdict saying
that the legacy Weinberg-angle role is replaced by an explicit strict
successor semantics?
```

## Result

The route remains negative on the current repo state:

```text
CURRENT_STRICT_SIDE_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_REPLACED_BRANCH_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P71
```

## What was checked

`P71` checks only the replaced branch. It asks whether the current strict-side
source set exports a statement that jointly contains:

1. the legacy Weinberg role or old formula,
2. a strict-side successor object or strict-side observable chain,
3. an explicit replacement/supersession verdict.

## Why it fails

On the current repo state:

1. the strict side does export `sin2_theta_w_mz` as a real candidate object,
2. the retained branch is already closed negatively,
3. but none of the current strict-side sources exports an explicit statement
   that the old legacy Weinberg role is replaced by a named strict successor
   semantics,
4. therefore the replaced branch remains non-discharged.

## Real reduction after `P71`

`P71` does not prove that the replaced branch is impossible forever.

It proves only the stronger current-repo-state result:

```text
the current strict-side source set exports no explicit replaced-branch verdict
for the legacy Weinberg-angle role
```

So the claim-specific Weinberg frontier is no longer:

```text
retained branch vs replaced branch
```

It is now:

```text
replaced branch only
```

still missing:

```text
explicit_strict_side_replaced_verdict_for_the_legacy_weinberg_angle_role_by_an_explicit_strict_successor_semantics
```

## What P71 does not claim

`P71` does not claim:

- that a replaced verdict can never exist,
- that the retained branch is globally impossible forever,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that the replaced branch is still absent on the
   current repo state,
2. then decide whether to refine the replaced branch into narrower successor
   sub-branches,
3. while keeping retained branch closed negatively on the current repo state.
