# F16 Legacy Fine-Structure Semantic-Transfer Refinement Packet

Status: `F16_EXECUTED_LEGACY_FINE_STRUCTURE_SEMANTIC_TRANSFER_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F15/P82/N87`, the retained-side frontier for the legacy fine-structure
role is narrowed to one missing object:

```text
explicit_legacy_to_strict_semantic_transfer_verdict_identifying_alpha_em_inv_mz
as_the_retained_strict_side_successor_of_the_legacy_fine_structure_role
```

`F16` asks the next honest refinement question:

```text
does that remaining semantic-transfer blocker split into
1. an explicit textual retained-successor verdict,
and
2. an explicit object-lineage-upgrade verdict that elevates the existing
   alpha_em_inv_mz candidate chain into actual retained semantic transfer?
```

## Result

`F16` establishes the following refinement on the current repo state:

1. the strict side already exports a real fine-structure candidate object
   `alpha_em_inv_mz`,
2. the repo also exports a real strict-side candidate chain around that object:
   `QW-2068` registry -> `QW-2069` derived package -> `QW-2098` strict gate,
3. `QW-2094` exports consistency markers for that candidate chain,
4. but that candidate chain is still weaker than an explicit retained
   semantic-transfer verdict,
5. therefore the remaining semantic-transfer blocker is now narrowed to two
   still-missing sub-blockers:
   an explicit textual retained-successor verdict and an explicit object-lineage
   upgrade verdict.

## Why this follows

The split is forced by current repo evidence:

1. `P82/N87` already show that a strict-side candidate object exists while a
   semantic-transfer verdict does not,
2. the repo exports a real object-side chain for `alpha_em_inv_mz`,
3. `QW-2094` confirms consistency markers for that chain,
4. but none of the current strict-side sources upgrades that chain into an
   explicit retained-role verdict,
5. so the next honest refinement is:
   textual retained-successor verdict vs object-lineage-upgrade verdict.

## Real reduction after `F16`

The retained semantic-transfer frontier is no longer:

```text
one generic semantic-transfer blocker
```

It is now:

```text
two narrower semantic-transfer sub-blockers
```

namely:

1. `explicit_textual_retained_successor_verdict_identifying_alpha_em_inv_mz_as_the_retained_strict_side_successor_of_the_legacy_fine_structure_role`
2. `explicit_object_lineage_upgrade_verdict_elevating_the_existing_alpha_em_inv_mz_candidate_chain_into_retained_strict_side_fine_structure_role_transfer`

## What F16 does not claim

`F16` does not claim:

- that either sub-blocker is already discharged,
- that the replaced branch is already discharged,
- that the `alpha_em_inv_mz` candidate chain itself already settles retained
  semantics,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. probe directly whether the current repo exports either of those two
   semantic-transfer sub-branches,
2. while keeping the replaced branch separate,
3. and without silently promoting the `alpha_em_inv_mz` candidate chain into
   full retained semantic transfer.
