# P64 Legacy Weinberg Role Strict-Side Branch Probe

Status: `P64_EXECUTED_LEGACY_WEINBERG_ROLE_STRICT_SIDE_BRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F6`, the missing strict-side verdict for the legacy Weinberg-angle role
is reduced to two narrower branches:

1. retained branch,
2. replaced branch.

`P64` asks the next honest question:

```text
does the current repo already export either branch for the legacy
Weinberg-angle role?
```

## Result

The route is negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_NEITHER_RETAINED_NOR_REPLACED_STRICT_SIDE_BRANCH_FOR_THE_LEGACY_WEINBERG_ANGLE_ROLE_AFTER_P64
```

## Why it fails

`P64` confirms all of the following at once:

1. the retained branch is not exported,
2. the replaced branch is not exported,
3. therefore the repo still has no claim-specific strict-side verdict for the
   legacy Weinberg-angle role.

## Real reduction after `P64`

`P64` does not say that one of those branches is impossible forever.

It says only this stronger current-repo-state result:

```text
the repo exports neither a retained verdict nor a replaced verdict for the
legacy Weinberg-angle role
```

So the frontier is no longer:

```text
one missing strict-side Weinberg verdict
```

It is now:

```text
two explicit missing branch verdicts
```

## What P64 does not claim

`P64` does not claim:

- theorem-level proof that the legacy Weinberg-angle role is retained,
- theorem-level proof that the legacy Weinberg-angle role is replaced,
- theorem-level proof that no future branch verdict can ever exist,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. attack one of the two branches,
2. most naturally the retained branch first,
3. while keeping non-transfer explicit.
