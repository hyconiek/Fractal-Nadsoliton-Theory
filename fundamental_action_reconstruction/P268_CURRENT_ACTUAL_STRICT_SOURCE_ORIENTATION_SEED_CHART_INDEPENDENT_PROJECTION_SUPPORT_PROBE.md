# P268 Current Actual Strict Source Orientation Seed Chart-Independent Projection Support Probe

Status: `P268_EXECUTED_CURRENT_ACTUAL_STRICT_SOURCE_ORIENTATION_SEED_CHART_INDEPENDENT_PROJECTION_SUPPORT_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P268` tests whether the current repo really exports the chart-independent
source-seed projection support packet introduced in `F177`, while keeping the
result:

1. below discharge of `strict_source_orientation_seed_target_v1`,
2. below actual chart-independent seed object export,
3. below actual pair-indexed population anchor,
4. below actual theta values,
5. below actual `E_orient`,
6. below admissible `S_sel_int`,
7. below strict-core selector closure,
8. below ToE closure.

## What P268 checks

`P268` checks only:

1. the first `T26` component remains future-only and undischarged,
2. the one-sided local descent support packet from `N287` remains exported,
3. the actual basis-free source-topology support layer from `N256/N257/N258`
   remains exported,
4. the local support packet can therefore be honestly paired with one
   chart-independent projection support route into that basis-free layer,
5. no actual chart-independent seed object, pair-indexing, theta values,
   `E_orient`, admissible `S_sel_int`, or closure claim is made.

## Result

`P268` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_SOURCE_ORIENTATION_SEED_CHART_INDEPENDENT_PROJECTION_SUPPORT_PACKET_BELOW_COMPONENT_DISCHARGE_AFTER_P268
```

This means:

1. the first `T26` component now has one stronger support packet than `N287`,
2. the support is stronger because it reaches the already actual basis-free
   source-topology layer,
3. the component still remains undischarged.

## Hard limits

`P268` does not establish:

1. discharge of `strict_source_orientation_seed_target_v1`,
2. actual chart-independent seed object export,
3. actual pair-indexed population anchor,
4. actual theta values,
5. actual populated basis-pair instance,
6. actual `E_orient`,
7. admissible `S_sel_int`,
8. actual strict-core selector closure,
9. actual ToE closure.
