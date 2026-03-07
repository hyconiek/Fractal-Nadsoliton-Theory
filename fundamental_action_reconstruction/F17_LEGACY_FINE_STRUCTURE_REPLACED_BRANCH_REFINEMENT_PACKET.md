# F17 Legacy Fine-Structure Replaced-Branch Refinement Packet

Status: `F17_EXECUTED_LEGACY_FINE_STRUCTURE_REPLACED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N90`, the only remaining claim-specific frontier for the legacy
fine-structure role is:

```text
explicit_strict_side_replaced_verdict_for_the_legacy_fine_structure_role_by_an_explicit_strict_successor_semantics
```

`F17` asks the next honest refinement question:

```text
does that one replaced-branch blocker split into
1. an explicit object-successor verdict around alpha_em_inv_mz,
and
2. an explicit method-successor-semantics verdict around
   qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r?
```

## Result

`F17` establishes the following refinement on the current repo state:

1. the strict side already exports a real successor-candidate object
   `alpha_em_inv_mz`,
2. the strict side also exports a real successor-candidate method lineage
   `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r`,
3. but neither candidate is yet accompanied by an explicit replaced-branch
   verdict for the legacy fine-structure role,
4. therefore the one replaced-branch blocker is now narrowed to two
   still-missing sub-blockers:
   an object-successor verdict and a method-successor-semantics verdict.

## Why this follows

The split is forced by current repo evidence:

1. `N90` already shows that the retained branch is closed negatively and the
   fine-structure claim-specific frontier has passed to the replaced branch,
2. `QW-2069/QW-2098/Release 4.9` do export a real strict-side object
   `alpha_em_inv_mz`,
3. those same sources also export the named derivation method
   `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r`,
4. therefore the narrowest honest refinement is:
   object-successor verdict vs method-successor-semantics verdict.

## Real reduction after `F17`

The replaced-branch frontier is no longer:

```text
one generic replaced-branch verdict
```

It is now:

```text
two narrower replaced-branch sub-blockers
```

namely:

1. `explicit_object_successor_verdict_identifying_alpha_em_inv_mz_as_the_strict_side_successor_object_replacing_the_legacy_fine_structure_role`
2. `explicit_method_successor_semantics_verdict_identifying_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_as_the_strict_side_successor_semantics_replacing_the_legacy_fine_structure_role`

## What F17 does not claim

`F17` does not claim:

- that either replaced sub-branch is already discharged,
- that `alpha_em_inv_mz` automatically replaces the legacy fine-structure role,
- that the `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r` chain
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
