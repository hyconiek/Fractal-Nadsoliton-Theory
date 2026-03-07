# P91 Strict-Side Textual Method Successor Probe For Legacy Fine-Structure Role

Status: `P91_EXECUTED_STRICT_SIDE_TEXTUAL_METHOD_SUCCESSOR_PROBE_FOR_LEGACY_FINE_STRUCTURE_ROLE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F19/P90/N95`, the next direct method-side question is:

```text
does the current strict-side method source set export an explicit textual
verdict identifying qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r as the
strict-side successor semantics replacing the legacy fine-structure role?
```

## Result

The route remains negative on the current repo state:

```text
CURRENT_STRICT_SIDE_METHOD_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_TEXTUAL_METHOD_SUCCESSOR_VERDICT_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P91
```

## What was checked

`P91` checks only the textual method-successor sub-branch. It asks whether the
current strict-side method source set exports a statement that jointly contains:

1. the legacy fine-structure role or old formula,
2. the strict-side method `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r`,
3. an explicit textual replacement/method-successor verdict.

## Why it fails

On the current repo state:

1. the method chain around `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r` is
   real,
2. the method-successor branch remains open after `P90`,
3. but no current strict-side method source exports an explicit textual verdict
   saying that `qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r` is the
   successor semantics replacing the legacy fine-structure role,
4. therefore the textual method-successor path remains non-discharged.

## Real reduction after `P91`

The method-successor branch is no longer blocked by two equally wide method-side
sub-branches.

It is now narrowed further:

```text
the textual method-successor sub-branch is negative on the current repo state
```

so the remaining method-side blocker is:

```text
explicit_method_lineage_upgrade_verdict_elevating_the_existing_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_chain_into_replacement_semantics_for_the_legacy_fine_structure_role
```

## What P91 does not claim

`P91` does not claim:

- that the method-lineage-upgrade sub-branch is impossible forever,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that the textual method-successor path is absent on
   the current repo state,
2. then attack the remaining method-lineage-upgrade sub-branch,
3. without reopening retained-side or object-side semantics.
