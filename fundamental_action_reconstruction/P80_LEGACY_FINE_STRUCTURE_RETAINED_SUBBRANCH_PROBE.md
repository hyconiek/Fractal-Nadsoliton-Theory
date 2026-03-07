# P80 Legacy Fine-Structure Retained Subbranch Probe

Status: `P80_EXECUTED_LEGACY_FINE_STRUCTURE_RETAINED_SUBBRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F14`, the strongest direct retained-side question is:

```text
does the current repo already export either literal retention or
role-equivalence retention for the legacy fine-structure role?
```

## Result

The route is negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_NEITHER_LITERAL_NOR_ROLE_EQUIVALENCE_RETAINED_SUBBRANCH_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P80
```

## Why it fails

`P80` confirms all of the following:

1. `F14` already reduces the retained branch to literal-retention and
   role-equivalence-retention subbranches,
2. the current repo exports neither literal strict-side retention of the old
   fine-structure formula,
3. nor explicit strict-side role-equivalence retention for the legacy
   fine-structure role.

## Real reduction after `P80`

The retained branch is no longer:

```text
one missing retained verdict
```

It is now:

```text
two retained-subbranch blockers
```

namely:

1. `explicit_strict_side_literal_retention_of_alpha_em_inverse_equals_alpha_geo_over_2beta_tors_times_1_minus_beta_tors`
2. `explicit_strict_side_role_equivalence_verdict_for_the_legacy_fine_structure_role`

## What P80 does not claim

`P80` does not claim:

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
