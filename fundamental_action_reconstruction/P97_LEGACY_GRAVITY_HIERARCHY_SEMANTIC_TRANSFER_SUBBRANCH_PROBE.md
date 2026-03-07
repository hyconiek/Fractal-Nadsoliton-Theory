# P97 Legacy Gravity-Hierarchy Semantic-Transfer Subbranch Probe

Status: `P97_EXECUTED_LEGACY_GRAVITY_HIERARCHY_SEMANTIC_TRANSFER_SUBBRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F23`, the next honest direct question is:

```text
does the current repo already export either of the two narrowed retained-side
semantic-transfer sub-branches for the legacy gravity-hierarchy role?
```

## Result

The route remains negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_NEITHER_TEXTUAL_SUCCESSOR_NOR_OBJECT_LINEAGE_UPGRADE_TRANSFER_VERDICT_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P97
```

## What was checked

`P97` checks two different narrowed sub-branches and keeps them separate:

1. whether the strict-side materials export an explicit textual verdict saying
   that `gravity_hierarchy_beta20` is the retained successor of the old legacy
   gravity-hierarchy role,
2. whether the repo exports an explicit lineage-upgrade verdict turning the
   existing `gravity_hierarchy_beta20` candidate chain into actual retained
   semantic transfer.

## Why it fails

On the current repo state:

1. the strict-side candidate object exists,
2. the strict-side candidate chain exists,
3. but neither a textual retained-successor verdict nor an explicit
   object-lineage-upgrade verdict is exported,
4. therefore the retained-side semantic-transfer route remains non-discharged.

## Real reduction after `P97`

The retained semantic-transfer frontier is no longer:

```text
one generic retained semantic-transfer blocker
```

It is now:

```text
two explicit still-missing sub-branch blockers
```

namely:

1. `explicit_textual_retained_successor_verdict_identifying_gravity_hierarchy_beta20_as_the_retained_strict_side_successor_of_the_legacy_gravity_hierarchy_role`
2. `explicit_object_lineage_upgrade_verdict_elevating_the_existing_gravity_hierarchy_beta20_candidate_chain_into_retained_strict_side_gravity_hierarchy_role_transfer`

## What P97 does not claim

`P97` does not claim:

- that either sub-branch is impossible forever,
- that the replaced branch is solved,
- that the `gravity_hierarchy_beta20` candidate chain is irrelevant,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that both semantic-transfer sub-branches remain
   absent on the current repo state,
2. then choose one of them for the next direct attack,
3. while keeping the replaced branch separate.
