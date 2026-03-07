# P94 Legacy Gravity-Hierarchy Retained Subbranch Probe

Status: `P94_EXECUTED_LEGACY_GRAVITY_HIERARCHY_RETAINED_SUBBRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F21`, the strongest direct retained-side question is:

```text
does the current repo already export either literal retention or
role-equivalence retention for the legacy gravity-hierarchy role?
```

## Result

The route is negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_NEITHER_LITERAL_NOR_ROLE_EQUIVALENCE_RETAINED_SUBBRANCH_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P94
```

## Why it fails

`P94` confirms all of the following:

1. `F21` already reduces the retained branch to literal-retention and
   role-equivalence-retention subbranches,
2. the current repo exports neither literal strict-side retention of the old
   gravity-hierarchy claim from `beta^N` scaling,
3. nor explicit strict-side role-equivalence retention for the legacy
   gravity-hierarchy role.

## Real reduction after `P94`

The retained branch is no longer:

```text
one missing retained verdict
```

It is now:

```text
two retained-subbranch blockers
```

namely:

1. `explicit_strict_side_literal_retention_of_exact_gravity_hierarchy_from_beta_to_the_N_scaling`
2. `explicit_strict_side_role_equivalence_verdict_for_the_legacy_gravity_hierarchy_role`

## What P94 does not claim

`P94` does not claim:

- theorem-level proof that literal retention is impossible forever,
- theorem-level proof that role-equivalence retention is impossible forever,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that both retained sub-branches remain absent,
2. then attack one retained sub-branch directly,
3. most naturally literal retention first.
