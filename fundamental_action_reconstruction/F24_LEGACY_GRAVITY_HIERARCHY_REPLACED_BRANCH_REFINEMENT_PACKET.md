# F24 Legacy Gravity-Hierarchy Replaced-Branch Refinement Packet

Status: `F24_EXECUTED_LEGACY_GRAVITY_HIERARCHY_REPLACED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N106`, the only remaining claim-specific frontier for the legacy
gravity-hierarchy role is:

```text
explicit_strict_side_replaced_verdict_for_the_legacy_gravity_hierarchy_role_by_an_explicit_strict_successor_semantics
```

`F24` asks the next honest refinement question:

```text
does that one replaced-branch blocker split into
1. an explicit object-successor verdict around gravity_hierarchy_beta20,
and
2. an explicit method-successor-semantics verdict around
   qw2115_micro_supported_beta_hierarchy_bridge?
```

## Result

`F24` establishes the following refinement on the current repo state:

1. the strict side already exports a real successor-candidate object
   `gravity_hierarchy_beta20`,
2. the strict side also exports a real successor-candidate method lineage
   `qw2115_micro_supported_beta_hierarchy_bridge`,
3. but neither candidate is yet accompanied by an explicit replaced-branch
   verdict for the legacy gravity-hierarchy role,
4. therefore the one replaced-branch blocker is now narrowed to two
   still-missing sub-blockers:
   an object-successor verdict and a method-successor-semantics verdict.

## Why this follows

The split is forced by current repo evidence:

1. `N106` already shows that the retained branch is closed negatively and the
   gravity-hierarchy claim-specific frontier has passed to the replaced
   branch,
2. `QW-2068/QW-2069/QW-2115/Release 4.9` do export a real strict-side object
   `gravity_hierarchy_beta20`,
3. `QW-2069/QW-2115` also export the named derivation method
   `qw2115_micro_supported_beta_hierarchy_bridge`,
4. therefore the narrowest honest refinement is:
   object-successor verdict vs method-successor-semantics verdict.

## Real reduction after `F24`

The replaced-branch frontier is no longer:

```text
one generic replaced-branch verdict
```

It is now:

```text
two narrower replaced-branch sub-blockers
```

namely:

1. `explicit_object_successor_verdict_identifying_gravity_hierarchy_beta20_as_the_strict_side_successor_object_replacing_the_legacy_gravity_hierarchy_role`
2. `explicit_method_successor_semantics_verdict_identifying_qw2115_micro_supported_beta_hierarchy_bridge_as_the_strict_side_successor_semantics_replacing_the_legacy_gravity_hierarchy_role`

## What F24 does not claim

`F24` does not claim:

- that either replaced sub-branch is already discharged,
- that `gravity_hierarchy_beta20` automatically replaces the legacy
  gravity-hierarchy role,
- that the `qw2115_micro_supported_beta_hierarchy_bridge` chain
  automatically upgrades into replacement semantics,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. probe directly whether the current repo exports either of those two
   replaced-branch sub-branches,
2. while keeping the retained branch closed negatively,
3. and without silently promoting strict candidate presence into replacement
   semantics.
