# P66 Strict-Side Literal Retention Probe For Legacy Weinberg Formula

Status: `P66_EXECUTED_STRICT_SIDE_LITERAL_RETENTION_PROBE_FOR_LEGACY_WEINBERG_FORMULA_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F7/P65/N68`, the next honest direct question is:

```text
do the current strict-side authoritative sources actually export the old
Weinberg formula or an algebraically identical literal form?
```

The strict-side authoritative source set used here is:

1. `RELEASE_4_9_TEXTBOOK_EN_PL.md`
2. `QW-2000`
3. `QW-2001`
4. `QW-2002`
5. `QW-2003`
6. `QW-2049`
7. `QW-2064`

## Result

The route is negative on the current repo state:

```text
CURRENT_STRICT_SIDE_AUTHORITATIVE_SOURCES_DO_NOT_EXPORT_LITERAL_RETENTION_OF_THE_LEGACY_WEINBERG_FORMULA_AFTER_P66
```

## What was checked

`P66` checks for literal or algebraically identical forms of the old legacy
Weinberg formula, including:

1. `sin^2(theta_W)=alpha_geo/12`
2. `sin^2(theta_W)=4 ln 2 / 12`
3. `sin^2(theta_W)=ln 2 / 3`

and the corresponding spacing/TeX variants.

## Why it fails

On the current repo state, none of the strict-side authoritative sources
exports any of those literal or algebraically identical forms.

So `P66` does **not** say that the old Weinberg role is impossible on the
strict side forever.

It says only this stronger current-repo-state result:

```text
the current strict-side authoritative source set does not literally retain the
old legacy Weinberg formula
```

## Real reduction after `P66`

This closes the `literal retention` sub-branch negatively on the current repo
state.

So the retained frontier for the legacy Weinberg role is no longer:

```text
literal retention vs role-equivalence retention
```

It is now:

```text
role-equivalence retention only
```

## What P66 does not claim

`P66` does not claim:

- that role-equivalence retention is absent,
- that the retained branch is globally impossible forever,
- that the replaced branch is solved,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. attack the remaining retained sub-branch directly:
   `explicit_strict_side_role_equivalence_verdict_for_the_legacy_weinberg_angle_role`,
2. while keeping the replaced branch separate.
