# P93 Legacy Gravity-Hierarchy Role Strict-Side Branch Probe

Status: `P93_EXECUTED_LEGACY_GRAVITY_HIERARCHY_ROLE_STRICT_SIDE_BRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F20`, the strongest direct branch-level question is:

```text
does the current repo already export either the retained or the replaced
strict-side branch for the legacy gravity-hierarchy role?
```

## Result

The route is negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_NEITHER_RETAINED_NOR_REPLACED_STRICT_SIDE_BRANCH_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P93
```

## Why it fails

`P93` confirms all of the following:

1. `F20` already reduces the missing gravity-hierarchy verdict to retained vs
   replaced branches,
2. the current repo still exports no explicit strict-side retained verdict for
   the legacy gravity-hierarchy role,
3. and no explicit strict-side replaced verdict for the legacy
   gravity-hierarchy role.

## Real reduction after `P93`

The frontier is no longer:

```text
one claim-specific strict-side verdict for the legacy gravity-hierarchy role
```

It is now:

```text
two branch-specific blockers
```

namely:

1. `explicit_strict_side_retained_verdict_for_the_legacy_gravity_hierarchy_role`
2. `explicit_strict_side_replaced_verdict_for_the_legacy_gravity_hierarchy_role_by_an_explicit_strict_successor_semantics`

## What P93 does not claim

`P93` does not claim:

- theorem-level proof that the legacy gravity-hierarchy role is retained,
- theorem-level proof that it is replaced,
- theorem-level proof that no future claim-specific partition can ever exist,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that both strict-side branch options remain absent,
2. then attack one branch directly,
3. most naturally the retained branch first.
