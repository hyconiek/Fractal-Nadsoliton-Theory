# F19 Legacy Fine-Structure Method-Successor Semantics Refinement Packet

Status: `F19_EXECUTED_LEGACY_FINE_STRUCTURE_METHOD_SUCCESSOR_SEMANTICS_REFINEMENT_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `N94`, the only remaining replaced-side frontier for the legacy
fine-structure role is:

```text
explicit_method_successor_semantics_verdict_identifying_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_as_the_strict_side_successor_semantics_replacing_the_legacy_fine_structure_role
```

`F19` asks the next honest refinement question:

```text
does that one method-side blocker split into
1. an explicit textual method-successor-semantics verdict,
and
2. an explicit method-lineage-upgrade verdict elevating the existing
   qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r method chain into replacement
   semantics?
```

## Result

`F19` establishes the following refinement on the current repo state:

1. the strict side already exports a real named method chain
   `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r`,
2. that method chain is visible across `QW-2069`, `QW-2098`, and the
   `QW-2094` method-consistency check for `alpha_em_inv_mz`,
3. but that chain is still weaker than an explicit replacement-semantics
   verdict for the legacy fine-structure role,
4. therefore the one method-side blocker is now narrowed to two still-missing
   sub-blockers:
   an explicit textual method-successor-semantics verdict and an explicit
   method-lineage-upgrade verdict.

## Why this follows

The split is forced by current repo evidence:

1. `N94` already closes the object-successor branch negatively,
2. the remaining replaced-side frontier is method-only,
3. the repo does export a real named strict method chain,
4. but none of the current sources upgrades that chain into explicit
   replacement semantics for the legacy fine-structure role,
5. so the next honest refinement is:
   textual method-successor verdict vs method-lineage-upgrade verdict.

## Real reduction after `F19`

The method-successor frontier is no longer:

```text
one generic method-successor-semantics verdict
```

It is now:

```text
two narrower method-successor sub-blockers
```

namely:

1. `explicit_textual_method_successor_semantics_verdict_identifying_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_as_the_strict_side_successor_semantics_replacing_the_legacy_fine_structure_role`
2. `explicit_method_lineage_upgrade_verdict_elevating_the_existing_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_chain_into_replacement_semantics_for_the_legacy_fine_structure_role`

## What F19 does not claim

`F19` does not claim:

- that either method-side sub-blocker is already discharged,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. probe directly whether the current repo exports either of those two
   method-side sub-branches,
2. without reopening retained-side or object-side semantics,
3. and without silently promoting method presence into replacement semantics.
