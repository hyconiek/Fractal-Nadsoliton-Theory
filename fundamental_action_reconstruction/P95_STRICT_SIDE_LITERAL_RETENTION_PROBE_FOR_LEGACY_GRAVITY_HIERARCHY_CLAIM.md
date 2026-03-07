# P95 Strict-Side Literal Retention Probe For Legacy Gravity-Hierarchy Claim

Status: `P95_EXECUTED_STRICT_SIDE_LITERAL_RETENTION_PROBE_FOR_LEGACY_GRAVITY_HIERARCHY_CLAIM_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F21/P94/N101`, the next honest direct question is:

```text
do the current strict-side authoritative sources actually export the old
gravity-hierarchy claim or an algebraically identical literal form?
```

The strict-side authoritative source set used here is:

1. `RELEASE_4_9_TEXTBOOK_EN_PL.md`
2. `QW-2069`
3. `QW-2094`
4. `QW-2115`
5. `A8`

## Result

The route is negative on the current repo state:

```text
CURRENT_STRICT_SIDE_AUTHORITATIVE_SOURCES_DO_NOT_EXPORT_LITERAL_RETENTION_OF_THE_LEGACY_GRAVITY_HIERARCHY_CLAIM_AFTER_P95
```

## What was checked

`P95` checks for literal or algebraically identical forms of the old legacy
gravity-hierarchy claim, including:

1. `exact gravity hierarchy from beta^N scaling`
2. `G_eff(N)=G_0*beta_tors^N`
3. `G(N=20)/G(N=0)=beta_tors^20=10^-40`

and the corresponding spacing/TeX variants.

## Why it fails

On the current repo state, none of the strict-side authoritative sources
exports any of those literal or algebraically identical forms.

So `P95` does **not** say that the old gravity-hierarchy role is impossible on
the strict side forever.

It says only this stronger current-repo-state result:

```text
the current strict-side authoritative source set does not literally retain the
old legacy gravity-hierarchy claim
```

## Real reduction after `P95`

This closes the `literal retention` sub-branch negatively on the current repo
state.

So the retained frontier for the legacy gravity-hierarchy role is no longer:

```text
literal retention vs role-equivalence retention
```

It is now:

```text
role-equivalence retention only
```

## What P95 does not claim

`P95` does not claim:

- that role-equivalence retention is absent,
- that the retained branch is globally impossible forever,
- that the replaced branch is solved,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. attack the remaining retained-side sub-branch directly:
   `explicit_strict_side_role_equivalence_verdict_for_the_legacy_gravity_hierarchy_role`,
2. while keeping the replaced branch separate.
