# P267 Current Actual Strict Source Orientation Seed One-Sided Descent Support Probe

Status: `P267_EXECUTED_CURRENT_ACTUAL_STRICT_SOURCE_ORIENTATION_SEED_ONE_SIDED_DESCENT_SUPPORT_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P267` tests whether the current repo really exports the stronger one-sided
local descent support packet introduced in `F176`, while keeping the result:

1. below discharge of `strict_source_orientation_seed_target_v1`,
2. below actual pair-indexed population anchor,
3. below actual theta values,
4. below actual `E_orient`,
5. below admissible `S_sel_int`,
6. below strict-core selector closure,
7. below ToE closure.

## What P267 checks

`P267` checks only:

1. the first `T26` component remains future-only and undischarged,
2. the stronger local branch-polarity support packet from `N286` remains
   exported,
3. the local derivative formula from `F163` is continuous for `d >= 0`,
4. the strict negativity at `d = 0` can therefore be extended to one positive
   local forward interval,
5. the protected positive source branch from `N249` remains exported,
6. these two facts can therefore be honestly combined into one forward
   half-branch on which `K_strict_gate(d)` stays positive and strictly
   decreases away from the source origin,
7. no chart-independent seed object, pair-indexing, actual theta values,
   actual `E_orient`, admissible `S_sel_int`, or closure claim is made.

## Result

`P267` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_SOURCE_ORIENTATION_SEED_ONE_SIDED_DESCENT_SUPPORT_PACKET_BELOW_COMPONENT_DISCHARGE_AFTER_P267
```

This means:

1. the first `T26` component now has one stronger support packet than `N286`,
2. the support is stronger because it reaches one distinguished positive
   forward half-branch with strict descent,
3. the component still remains undischarged.

## Hard limits

`P267` does not establish:

1. discharge of `strict_source_orientation_seed_target_v1`,
2. actual pair-indexed population anchor,
3. actual theta values,
4. actual populated basis-pair instance,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. actual strict-core selector closure,
8. actual ToE closure.
