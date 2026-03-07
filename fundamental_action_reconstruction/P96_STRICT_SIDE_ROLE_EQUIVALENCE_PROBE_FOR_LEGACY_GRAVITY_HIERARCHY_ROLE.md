# P96 Strict-Side Role-Equivalence Probe For Legacy Gravity-Hierarchy Role

Status: `P96_EXECUTED_STRICT_SIDE_ROLE_EQUIVALENCE_PROBE_FOR_LEGACY_GRAVITY_HIERARCHY_ROLE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F22`, the next honest direct question is:

```text
does the current repo already export an explicit legacy-to-strict
role-equivalence verdict for the legacy gravity-hierarchy role?
```

## Result

The route is mixed but still negative at the theorem-relevant level:

```text
CURRENT_REPO_EXPORTS_STRICT_SIDE_GRAVITY_HIERARCHY_BETA20_CANDIDATE_BUT_NO_EXPLICIT_LEGACY_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_VERDICT_AFTER_P96
```

## What was checked

`P96` checks two different things and keeps them separate:

1. whether the strict side exports a real candidate object carrying the modern
   gravity-hierarchy observable label,
2. whether the repo also exports an explicit semantic-transfer verdict saying
   that this object is the retained role-equivalent successor of the old
   legacy gravity-hierarchy semantics.

## Why it only partially succeeds

On the current repo state:

1. `QW-2068`, `QW-2069`, and `QW-2115` together do export a real strict-side
   candidate object `gravity_hierarchy_beta20`,
2. but none of the current strict-side sources exports an explicit verdict
   that `gravity_hierarchy_beta20` is the retained role-equivalent successor
   of the legacy gravity-hierarchy semantics,
3. therefore the retained-side role-equivalence branch is narrowed but not yet
   discharged.

## Real reduction after `P96`

The retained frontier for the legacy gravity-hierarchy role is no longer:

```text
generic role-equivalence retention blocker
```

It is now:

```text
one explicit semantic-transfer blocker attached to gravity_hierarchy_beta20
```

namely:

```text
explicit_legacy_to_strict_semantic_transfer_verdict_identifying_gravity_hierarchy_beta20_as_the_retained_strict_side_successor_of_the_legacy_gravity_hierarchy_role
```

## What P96 does not claim

`P96` does not claim:

- that the retained branch is already discharged,
- that the replaced branch is already discharged,
- that `gravity_hierarchy_beta20` automatically inherits all legacy
  semantics,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that the current repo exports a strict-side
   gravity-hierarchy candidate object but no explicit legacy-to-strict
   role-equivalence verdict,
2. keep the replaced branch separate,
3. do not silently treat `gravity_hierarchy_beta20` as if it already carried
   the full legacy `beta^N` semantics.
