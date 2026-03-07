# F10 Legacy Weinberg Replaced-Branch Refinement Packet

Status: `F10_EXECUTED_LEGACY_WEINBERG_REPLACED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P71/N74`, the only remaining claim-specific frontier for the legacy
Weinberg-angle role is:

```text
explicit_strict_side_replaced_verdict_for_the_legacy_weinberg_angle_role_by_an_explicit_strict_successor_semantics
```

`F10` asks the next honest refinement question:

```text
does that one replaced-branch blocker split into
1. an explicit object-successor verdict around sin2_theta_w_mz,
and
2. an explicit method-successor-semantics verdict around
   qw2098_sin2_from_nonanchor_ew_pole_chain?
```

## Result

`F10` establishes the following refinement on the current repo state:

1. the strict side already exports a real successor-candidate object
   `sin2_theta_w_mz`,
2. the strict side also exports a real successor-candidate method lineage
   `qw2098_sin2_from_nonanchor_ew_pole_chain`,
3. but neither candidate is yet accompanied by an explicit replaced-branch
   verdict for the legacy Weinberg-angle role,
4. therefore the one replaced-branch blocker is now narrowed to two still-missing
   sub-blockers:
   an object-successor verdict and a method-successor-semantics verdict.

## Why this follows

The split is forced by current repo evidence:

1. `P71/N74` already show that the replaced branch is still absent as a whole,
2. `QW-2069/QW-2098/Release 4.9` do export a real strict-side object
   `sin2_theta_w_mz`,
3. those same sources also export the named derivation method
   `qw2098_sin2_from_nonanchor_ew_pole_chain`,
4. therefore the narrowest honest refinement is:
   object-successor verdict vs method-successor-semantics verdict.

## Real reduction after `F10`

The replaced-branch frontier is no longer:

```text
one generic replaced-branch verdict
```

It is now:

```text
two narrower replaced-branch sub-blockers
```

namely:

1. `explicit_object_successor_verdict_identifying_sin2_theta_w_mz_as_the_strict_side_successor_object_replacing_the_legacy_weinberg_angle_role`
2. `explicit_method_successor_semantics_verdict_identifying_qw2098_sin2_from_nonanchor_ew_pole_chain_as_the_strict_side_successor_semantics_replacing_the_legacy_weinberg_angle_role`

## What F10 does not claim

`F10` does not claim:

- that either replaced sub-branch is already discharged,
- that `sin2_theta_w_mz` automatically replaces the legacy Weinberg role,
- that the `qw2098` chain automatically upgrades into replacement semantics,
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
