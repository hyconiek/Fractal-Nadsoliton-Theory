# F21 Legacy Gravity-Hierarchy Retained Branch Refinement Packet

Status: `F21_EXECUTED_LEGACY_GRAVITY_HIERARCHY_RETAINED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F20/P93/N100`, one narrow branch remains open:

```text
explicit_strict_side_retained_verdict_for_the_legacy_gravity_hierarchy_role
```

`F21` asks the next honest refinement question:

```text
what are the two narrowest ways the retained branch could be discharged on the
current repo state?
```

## Result

`F21` establishes the following refinement:

1. one possible retained-branch discharge would be an explicit strict-side
   literal retention of the old gravity-hierarchy claim
   `exact gravity hierarchy from beta^N scaling`,
2. the other possible retained-branch discharge would be an explicit
   strict-side role-equivalence verdict showing that some strict-side object
   carries the same gravity-hierarchy role even if the literal old `beta^N`
   claim is not retained,
3. the current repo exports neither one yet.

## Real reduction after F21

So the retained branch is no longer:

```text
one missing retained verdict
```

It is now:

```text
two narrower retained sub-branches
```

namely:

1. `explicit_strict_side_literal_retention_of_exact_gravity_hierarchy_from_beta_to_the_N_scaling`
2. `explicit_strict_side_role_equivalence_verdict_for_the_legacy_gravity_hierarchy_role`

## Why this follows

The split is forced by current repo evidence:

1. the old legacy gravity-hierarchy role is explicit in `QW-2005` and the
   legacy TeX package,
2. the strict-side materials do not currently export the old `beta^N` claim,
3. `P93/N100` already show that the repo has no overall retained verdict,
4. therefore the next honest decomposition is:
   literal retention vs role-equivalence retention.

## What F21 does not claim

`F21` does not claim:

- that the literal retention branch is present,
- that the role-equivalence branch is present,
- that either retained sub-branch is impossible forever,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. test whether the current repo exports either retained sub-branch,
2. or formalize theorem-level that both retained sub-branches remain absent.
