# F13 Legacy Fine-Structure Role Strict-Side Partition Refinement Packet

Status: `F13_EXECUTED_LEGACY_FINE_STRUCTURE_ROLE_STRICT_SIDE_PARTITION_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F5/P63/N66`, one claim-specific partition gap is now singled out:

```text
explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_fine_structure_role
```

`F13` asks the next honest refinement question:

```text
what are the two narrowest ways that one missing strict-side fine-structure
verdict could in principle be discharged on the current repo state?
```

## Result

`F13` establishes the following refinement:

1. one possible discharge would be an explicit strict-side `retained` verdict
   for the legacy fine-structure role,
2. the other possible discharge would be an explicit strict-side `replaced`
   verdict for the legacy fine-structure role by an explicit strict successor
   semantics,
3. the current repo exports neither one yet.

## Real reduction after F13

So the missing object is no longer:

```text
one claim-specific strict-side verdict for the legacy fine-structure role
```

It is now:

```text
two narrower verdict branches
```

namely:

1. `explicit_strict_side_retained_verdict_for_the_legacy_fine_structure_role`
2. `explicit_strict_side_replaced_verdict_for_the_legacy_fine_structure_role_by_an_explicit_strict_successor_semantics`

## Why this follows

The split is forced by the current repo state:

1. `P62/N65` already rule out silent transfer of the legacy fine-structure
   role onto `K_strict_gate`,
2. `P63/N66` already rule out any claim-specific strict-side verdict for the
   fine-structure role,
3. therefore the next honest decomposition is exactly:
   retained branch vs replaced branch.

## What F13 does not claim

`F13` does not claim:

- that the legacy fine-structure role is retained on the strict side,
- that the legacy fine-structure role is replaced on the strict side,
- that either branch is impossible forever,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. test whether the current repo already exports either the retained branch
   or the replaced branch for the legacy fine-structure role,
2. or formalize theorem-level that both narrow branches remain absent.
