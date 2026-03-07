# F14 Legacy Fine-Structure Retained Branch Refinement Packet

Status: `F14_EXECUTED_LEGACY_FINE_STRUCTURE_RETAINED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F13/P79/N84`, one narrow branch remains open:

```text
explicit_strict_side_retained_verdict_for_the_legacy_fine_structure_role
```

`F14` asks the next honest refinement question:

```text
what are the two narrowest ways the retained branch could be discharged on the
current repo state?
```

## Result

`F14` establishes the following refinement:

1. one possible retained-branch discharge would be an explicit strict-side
   literal retention of the old formula
   `alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)`,
2. the other possible retained-branch discharge would be an explicit
   strict-side role-equivalence verdict showing that some strict-side object
   carries the same fine-structure role even if the literal old formula is not
   retained,
3. the current repo exports neither one yet.

## Real reduction after F14

So the retained branch is no longer:

```text
one missing retained verdict
```

It is now:

```text
two narrower retained sub-branches
```

namely:

1. `explicit_strict_side_literal_retention_of_alpha_em_inverse_equals_alpha_geo_over_2beta_tors_times_1_minus_beta_tors`
2. `explicit_strict_side_role_equivalence_verdict_for_the_legacy_fine_structure_role`

## Why this follows

The split is forced by current repo evidence:

1. the old legacy fine-structure role is explicit in legacy sources,
2. the strict-side materials do not currently export the old formula,
3. `P79/N84` already show that the repo has no overall retained verdict,
4. therefore the next honest decomposition is:
   literal retention vs role-equivalence retention.

## What F14 does not claim

`F14` does not claim:

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
