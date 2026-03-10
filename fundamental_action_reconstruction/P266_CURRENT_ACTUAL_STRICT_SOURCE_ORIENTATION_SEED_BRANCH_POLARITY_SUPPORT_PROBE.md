# P266 Current Actual Strict Source Orientation Seed Branch Polarity Support Probe

Status: `P266_EXECUTED_CURRENT_ACTUAL_STRICT_SOURCE_ORIENTATION_SEED_BRANCH_POLARITY_SUPPORT_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P266` tests whether the current repo really exports the stronger local
branch-polarity support packet introduced in `F175`, while keeping the result:

1. below discharge of `strict_source_orientation_seed_target_v1`,
2. below actual pair-indexed population anchor,
3. below actual theta values,
4. below actual `E_orient`,
5. below admissible `S_sel_int`,
6. below strict-core selector closure,
7. below ToE closure.

## What P266 checks

`P266` checks only:

1. the first `T26` component remains future-only and undischarged,
2. the actual narrow support packet from `N285` remains exported,
3. the actual source-side barrier-protected sign witness from `N249` remains
   exported,
4. the local derivative datum from `F163/N273` remains exported and strictly
   negative at the source origin,
5. these two facts can therefore be honestly read together as one protected
   positive local branch plus one immediate descent polarity datum,
6. the actual source-topology selector support family from `N256/N257/N258`
   remains exported,
7. no chart-independent seed object, pair-indexing, actual theta values,
   actual `E_orient`, admissible `S_sel_int`, or closure claim is made.

## Result

`P266` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_SOURCE_ORIENTATION_SEED_BRANCH_POLARITY_SUPPORT_PACKET_BELOW_COMPONENT_DISCHARGE_AFTER_P266
```

This means:

1. the first `T26` component now has one stronger local support packet than
   `N285`,
2. the support is stronger because it combines protected branch positivity
   with strict local derivative polarity,
3. the component still remains undischarged.

## Hard limits

`P266` does not establish:

1. discharge of `strict_source_orientation_seed_target_v1`,
2. actual pair-indexed population anchor,
3. actual theta values,
4. actual populated basis-pair instance,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. actual strict-core selector closure,
8. actual ToE closure.
