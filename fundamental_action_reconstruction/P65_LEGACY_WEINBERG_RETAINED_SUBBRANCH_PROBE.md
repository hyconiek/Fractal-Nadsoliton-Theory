# P65 Legacy Weinberg Retained Subbranch Probe

Status: `P65_EXECUTED_LEGACY_WEINBERG_RETAINED_SUBBRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F7`, the retained branch for the legacy Weinberg-angle role is reduced
to two narrower sub-branches:

1. literal retention of `sin^2(theta_W)=alpha_geo/12`,
2. explicit strict-side role-equivalence retention.

`P65` asks the next honest question:

```text
does the current repo already export either retained sub-branch?
```

## Result

The route is negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_NEITHER_LITERAL_NOR_ROLE_EQUIVALENCE_RETAINED_SUBBRANCH_FOR_THE_LEGACY_WEINBERG_ANGLE_ROLE_AFTER_P65
```

## Why it fails

`P65` confirms all of the following at once:

1. no strict-side literal retention of
   `sin^2(theta_W)=alpha_geo/12` is currently exported,
2. no explicit strict-side role-equivalence retention verdict is currently
   exported,
3. therefore the retained branch for the legacy Weinberg-angle role remains
   non-discharged.

## Real reduction after `P65`

`P65` does not say that either retained sub-branch is impossible forever.

It says only this stronger current-repo-state result:

```text
the repo exports neither literal retention nor role-equivalence retention for
the legacy Weinberg-angle role
```

So the retained frontier is no longer:

```text
one missing retained branch
```

It is now:

```text
two explicit missing retained sub-branches
```

## What P65 does not claim

`P65` does not claim:

- theorem-level proof that literal retention can never exist,
- theorem-level proof that role-equivalence retention can never exist,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. attack one of the two retained sub-branches directly,
2. most naturally literal retention first,
3. while keeping the replaced branch separate.
