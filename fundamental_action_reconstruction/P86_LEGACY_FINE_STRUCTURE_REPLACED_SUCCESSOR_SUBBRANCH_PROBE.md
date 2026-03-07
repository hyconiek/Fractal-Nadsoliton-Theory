# P86 Legacy Fine-Structure Replaced Successor Sub-Branch Probe

Status: `P86_EXECUTED_LEGACY_FINE_STRUCTURE_REPLACED_SUCCESSOR_SUBBRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F17`, the replaced branch for the legacy fine-structure role has two
narrower sub-branches:

1. object-successor verdict,
2. method-successor-semantics verdict.

`P86` asks:

```text
does the current repo export either of those two replaced-side successor
sub-branches?
```

## Result

The route remains negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_STRICT_FINE_STRUCTURE_SUCCESSOR_CANDIDATES_BUT_NEITHER_OBJECT_NOR_METHOD_REPLACED_SUCCESSOR_VERDICT_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P86
```

## What was checked

`P86` checks only the replaced branch and only at the newly refined level:

1. whether any current strict-side source exports an explicit replaced
   object-successor verdict tying `alpha_em_inv_mz` to the legacy
   fine-structure role,
2. whether any current strict-side source exports an explicit replaced
   method-successor-semantics verdict tying
   `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r` to the legacy
   fine-structure role.

## Why it fails

On the current repo state:

1. the strict side does export `alpha_em_inv_mz`,
2. the strict side does export the named method
   `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r`,
3. but no current strict-side source exports an explicit replaced verdict on
   either sub-branch,
4. therefore the refined replaced branch remains non-discharged.

## Real reduction after `P86`

The replaced frontier is no longer:

```text
one generic replaced verdict
```

It is now:

```text
two narrower replaced successor sub-branches
```

still missing:

1. `explicit_object_successor_verdict_identifying_alpha_em_inv_mz_as_the_strict_side_successor_object_replacing_the_legacy_fine_structure_role`
2. `explicit_method_successor_semantics_verdict_identifying_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_as_the_strict_side_successor_semantics_replacing_the_legacy_fine_structure_role`

## What P86 does not claim

`P86` does not claim:

- that either sub-branch is impossible forever,
- that the retained branch reopens,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that both refined replaced sub-branches remain
   absent on the current repo state,
2. then attack one of them directly,
3. most naturally the object-successor verdict first.
