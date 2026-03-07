# P100 Legacy Gravity-Hierarchy Replaced-Successor Subbranch Probe

Status: `P100_EXECUTED_LEGACY_GRAVITY_HIERARCHY_REPLACED_SUCCESSOR_SUBBRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F24`, the next honest direct question is:

```text
does the current repo already export either of the two narrowed replaced-side
sub-branches for the legacy gravity-hierarchy role?
```

## Result

The route remains negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_STRICT_GRAVITY_HIERARCHY_SUCCESSOR_CANDIDATES_BUT_NEITHER_OBJECT_NOR_METHOD_REPLACED_SUCCESSOR_VERDICT_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P100
```

## What was checked

`P100` checks two different narrowed sub-branches and keeps them separate:

1. whether the strict-side materials export an explicit object-successor
   verdict saying that `gravity_hierarchy_beta20` replaces the old legacy
   gravity-hierarchy role,
2. whether the repo exports an explicit method-successor-semantics verdict
   turning `qw2115_micro_supported_beta_hierarchy_bridge` into actual
   replacement semantics for the same role.

## Why it fails

On the current repo state:

1. the strict-side candidate object exists,
2. the strict-side method chain exists,
3. but neither an object-successor replaced verdict nor an explicit
   method-successor-semantics verdict is exported,
4. therefore the replaced-side route remains non-discharged.

## Real reduction after `P100`

The replaced-side frontier is no longer:

```text
one generic replaced-branch blocker
```

It is now:

```text
two explicit still-missing sub-branch blockers
```

namely:

1. `explicit_object_successor_verdict_identifying_gravity_hierarchy_beta20_as_the_strict_side_successor_object_replacing_the_legacy_gravity_hierarchy_role`
2. `explicit_method_successor_semantics_verdict_identifying_qw2115_micro_supported_beta_hierarchy_bridge_as_the_strict_side_successor_semantics_replacing_the_legacy_gravity_hierarchy_role`

## What P100 does not claim

`P100` does not claim:

- that either sub-branch is impossible forever,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that both replaced-side sub-branches remain absent
   on the current repo state,
2. then choose one of them for the next direct attack,
3. while keeping retained and replaced branches separated.
