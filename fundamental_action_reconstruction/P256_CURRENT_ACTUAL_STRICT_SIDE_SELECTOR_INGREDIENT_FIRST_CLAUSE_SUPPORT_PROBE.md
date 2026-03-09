# P256 Current Actual Strict-Side Selector Ingredient First Clause Support Probe

Status: `P256_EXECUTED_CURRENT_ACTUAL_STRICT_SIDE_SELECTOR_INGREDIENT_FIRST_CLAUSE_SUPPORT_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P256` tests whether the current repo really exports the first-clause support
packet introduced in `F166`, while keeping the result:

1. below admissible `S_sel_int`,
2. below strict-core selector closure,
3. below ToE closure.

## What P256 checks

`P256` checks only:

1. the old genuine-source admission contract remains active,
2. the old first-clause obstruction remains active,
3. the repo now also exports actual source-topology support upstream of
   observer,
4. that support remains declared-scope only and kernel-split-safe,
5. no identity claim `tau_src_candidate_v1 == S_sel_int` is made,
6. no strict-core selector closure is claimed.

## Result

`P256` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_SIDE_SELECTOR_INGREDIENT_FIRST_CLAUSE_SUPPORT_PACKET_BELOW_ADMISSIBLE_S_SEL_INT_AFTER_P256
```

This means:

1. the strict-side frontier now contains one actual support packet for
   re-attacking the first admissibility clause,
2. the support packet is stronger than pure negative packaging,
3. the repo still does not justify saying that a genuine strict-core selector
   source already exists.

## Hard limits

`P256` does not establish:

1. admissible `S_sel_int`,
2. actual strict-core selector closure,
3. actual global `QW-2191` discharge,
4. actual ToE closure.
