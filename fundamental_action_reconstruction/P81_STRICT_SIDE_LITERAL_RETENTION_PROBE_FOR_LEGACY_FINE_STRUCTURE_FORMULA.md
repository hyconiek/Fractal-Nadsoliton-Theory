# P81 Strict-Side Literal Retention Probe For Legacy Fine-Structure Formula

Status: `P81_EXECUTED_STRICT_SIDE_LITERAL_RETENTION_PROBE_FOR_LEGACY_FINE_STRUCTURE_FORMULA_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F14/P80/N85`, the next honest direct question is:

```text
do the current strict-side authoritative sources actually export the old
fine-structure formula or an algebraically identical literal form?
```

The strict-side authoritative source set used here is:

1. `RELEASE_4_9_TEXTBOOK_EN_PL.md`
2. `QW-2069`
3. `QW-2094`
4. `QW-2098`

## Result

The route is negative on the current repo state:

```text
CURRENT_STRICT_SIDE_AUTHORITATIVE_SOURCES_DO_NOT_EXPORT_LITERAL_RETENTION_OF_THE_LEGACY_FINE_STRUCTURE_FORMULA_AFTER_P81
```

## What was checked

`P81` checks for literal or algebraically identical forms of the old legacy
fine-structure formula, including:

1. `alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)`
2. algebraically identical spacing/TeX variants of the same old formula.

## Why it fails

On the current repo state, none of the strict-side authoritative sources
exports any of those literal or algebraically identical forms.

So `P81` does **not** say that the old fine-structure role is impossible on the
strict side forever.

It says only this stronger current-repo-state result:

```text
the current strict-side authoritative source set does not literally retain the
old legacy fine-structure formula
```

## Real reduction after `P81`

This closes the `literal retention` sub-branch negatively on the current repo
state.

So the retained frontier for the legacy fine-structure role is no longer:

```text
literal retention vs role-equivalence retention
```

It is now:

```text
role-equivalence retention only
```

## What P81 does not claim

`P81` does not claim:

- that role-equivalence retention is absent,
- that the retained branch is globally impossible forever,
- that the replaced branch is solved,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. attack the remaining retained-side sub-branch directly:
   `explicit_strict_side_role_equivalence_verdict_for_the_legacy_fine_structure_role`,
2. while keeping the replaced branch separate.
