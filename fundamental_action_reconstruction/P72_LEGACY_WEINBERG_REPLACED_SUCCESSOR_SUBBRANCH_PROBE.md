# P72 Legacy Weinberg Replaced Successor Sub-Branch Probe

Status: `P72_EXECUTED_LEGACY_WEINBERG_REPLACED_SUCCESSOR_SUBBRANCH_PROBE_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `F10`, the replaced branch for the legacy Weinberg-angle role has two
narrower sub-branches:

1. object-successor verdict,
2. method-successor-semantics verdict.

`P72` asks:

```text
does the current repo export either of those two replaced-side successor
sub-branches?
```

## Result

The route remains negative on the current repo state:

```text
CURRENT_REPO_EXPORTS_STRICT_WEINBERG_SUCCESSOR_CANDIDATES_BUT_NEITHER_OBJECT_NOR_METHOD_REPLACED_SUCCESSOR_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P72
```

## What was checked

`P72` checks only the replaced branch and only at the newly refined level:

1. whether any current strict-side source exports an explicit replaced
   object-successor verdict tying `sin2_theta_w_mz` to the legacy
   Weinberg-angle role,
2. whether any current strict-side source exports an explicit replaced
   method-successor-semantics verdict tying
   `qw2098_sin2_from_nonanchor_ew_pole_chain` to the legacy Weinberg-angle
   role.

## Why it fails

On the current repo state:

1. the strict side does export `sin2_theta_w_mz`,
2. the strict side does export the named method
   `qw2098_sin2_from_nonanchor_ew_pole_chain`,
3. but no current strict-side source exports an explicit replaced verdict on
   either sub-branch,
4. therefore the refined replaced branch remains non-discharged.

## Real reduction after `P72`

The replaced frontier is no longer:

```text
one generic replaced verdict
```

It is now:

```text
two narrower replaced successor sub-branches
```

still missing:

1. `explicit_object_successor_verdict_identifying_sin2_theta_w_mz_as_the_strict_side_successor_object_replacing_the_legacy_weinberg_angle_role`
2. `explicit_method_successor_semantics_verdict_identifying_qw2098_sin2_from_nonanchor_ew_pole_chain_as_the_strict_side_successor_semantics_replacing_the_legacy_weinberg_angle_role`

## What P72 does not claim

`P72` does not claim:

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
