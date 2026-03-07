# F22 Legacy Gravity-Hierarchy Role-Equivalence Refinement Packet

Status: `F22_EXECUTED_LEGACY_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P95/N102`, the retained-side frontier for the legacy gravity-hierarchy
role is no longer literal retention. The only remaining retained-side
question is:

```text
explicit_strict_side_role_equivalence_verdict_for_the_legacy_gravity_hierarchy_role
```

`F22` asks the next honest refinement question:

```text
does that remaining role-equivalence frontier reduce to
1. the presence of some strict-side gravity-hierarchy candidate object,
or instead to
2. one still-missing explicit semantic-transfer verdict from legacy to strict?
```

## Result

`F22` establishes the following refinement on the current repo state:

1. the strict side already exports a real candidate object,
   `gravity_hierarchy_beta20`,
2. that candidate is not yet the same thing as an explicit retained-role
   verdict for the legacy gravity-hierarchy semantics,
3. therefore the remaining retained-side blocker is now narrowed to one
   explicit semantic-transfer / role-equivalence verdict linking the legacy
   role to the strict-side candidate object.

## Why this follows

The split is forced by current repo evidence:

1. `QW-2068` registers `gravity_hierarchy_beta20` as a target parameter,
2. `QW-2069` exports `gravity_hierarchy_beta20` as a strict-derived package
   entry,
3. `QW-2115` exports `gravity_hierarchy_beta20` as a strict internal gravity
   hierarchy bridge gate result,
4. but the current repo still exports no explicit verdict saying that this
   strict-side object is the retained role-equivalent successor of the old
   legacy gravity-hierarchy claim from `beta^N` scaling.

## Real reduction after `F22`

So the retained-side frontier is no longer:

```text
generic role-equivalence retention blocker
```

It is now:

```text
one explicit semantic-transfer blocker
```

namely:

```text
explicit_legacy_to_strict_semantic_transfer_verdict_identifying_gravity_hierarchy_beta20_as_the_retained_strict_side_successor_of_the_legacy_gravity_hierarchy_role
```

## What F22 does not claim

`F22` does not claim:

- that the retained branch is already discharged,
- that the replaced branch is already discharged,
- that `gravity_hierarchy_beta20` automatically inherits all legacy
  semantics,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. probe directly whether the current strict-side materials export that one
   missing semantic-transfer verdict,
2. while keeping the replaced branch separate and without silently
   transferring legacy semantics onto `K_strict_gate`.
