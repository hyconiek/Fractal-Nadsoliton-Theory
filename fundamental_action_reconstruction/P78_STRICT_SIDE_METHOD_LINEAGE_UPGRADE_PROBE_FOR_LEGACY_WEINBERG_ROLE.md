# P78 Strict-Side Method-Lineage Upgrade Probe For Legacy Weinberg Role

Status: `P78_EXECUTED_STRICT_SIDE_METHOD_LINEAGE_UPGRADE_PROBE_FOR_LEGACY_WEINBERG_ROLE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `P77/N80`, the last remaining method-side blocker for the legacy
Weinberg-angle role is:

```text
explicit_method_lineage_upgrade_verdict_elevating_the_existing
qw2098_sin2_from_nonanchor_ew_pole_chain_chain_into_replacement_semantics_for_the_legacy_weinberg_angle_role
```

`P78` asks:

```text
does the current strict-side method source set export an explicit
method-lineage-upgrade verdict elevating the qw2098 method chain into
replacement semantics for the legacy Weinberg-angle role?
```

## Result

The route remains negative on the current repo state:

```text
CURRENT_STRICT_SIDE_METHOD_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_METHOD_LINEAGE_UPGRADE_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P78
```

## What was checked

`P78` checks only the method-lineage-upgrade sub-branch. It asks whether the
current strict-side method source set exports a statement that jointly contains:

1. the legacy Weinberg-angle role or old formula,
2. the strict-side method chain `qw2098_sin2_from_nonanchor_ew_pole_chain`,
3. explicit lineage/consistency markers for that method chain,
4. an explicit upgrade into replacement semantics for the legacy
   Weinberg-angle role.

## Why it fails

On the current repo state:

1. the `qw2098_sin2_from_nonanchor_ew_pole_chain` method chain is real,
2. the textual method-successor sub-branch is already closed negatively by
   `P77/N80`,
3. but no current strict-side method source exports an explicit verdict
   upgrading the existing `qw2098` method chain into replacement semantics for
   the legacy Weinberg-angle role,
4. therefore the method-lineage-upgrade path remains non-discharged.

## Real reduction after `P78`

The method-side frontier is no longer:

```text
textual method-successor vs method-lineage upgrade
```

It is now:

```text
the full method-successor branch is negative on the current repo state
```

## What P78 does not claim

`P78` does not claim:

- that the method-lineage-upgrade path is impossible forever,
- that the full replaced branch is already theorem-level closed,
- selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is now:

1. formalize theorem-level that the full method-successor branch is closed
   negatively on the current repo state,
2. then combine that result with `N78` if one wants a full replaced-branch
   closure theorem,
3. without reopening retained-side or object-side semantics.
