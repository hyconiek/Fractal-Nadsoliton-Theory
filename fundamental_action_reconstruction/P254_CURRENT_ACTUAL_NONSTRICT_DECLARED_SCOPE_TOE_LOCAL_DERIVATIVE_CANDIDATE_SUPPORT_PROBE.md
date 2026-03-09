# P254 Current Actual Nonstrict Declared-Scope ToE Local Derivative Candidate Support Probe

Status: `P254_EXECUTED_CURRENT_ACTUAL_NONSTRICT_DECLARED_SCOPE_TOE_LOCAL_DERIVATIVE_CANDIDATE_SUPPORT_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P254` tests whether the current repo really exports the new candidate-support
packet introduced in `F164`, while keeping the result:

1. below selector closure,
2. below global `QW-2191` discharge,
3. below actual non-strict declared-scope ToE closure.

## What P254 checks

`P254` checks only:

1. the earlier non-strict declared-scope ToE preclosure support packet is
   still exported,
2. the future-only non-strict declared-scope ToE target is still exported,
3. the local derivative datum from `F163` is now integrated into one explicit
   candidate-support packet,
4. the derivative datum remains local, source-side only, and
   candidate-only,
5. no actual ToE closure is claimed from that packet.

## Result

`P254` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_NONSTRICT_DECLARED_SCOPE_TOE_LOCAL_DERIVATIVE_CANDIDATE_SUPPORT_PACKET_BELOW_CLOSURE_AFTER_P254
```

This means:

1. the non-strict declared-scope ToE lane now contains one additional
   packaged candidate ingredient,
2. that ingredient is still explicitly weaker than any closure discharge,
3. the current repo still does not justify actual ToE closure language.

## Hard limits

`P254` does not establish:

1. actual non-strict declared-scope ToE closure,
2. actual strict-core ToE closure,
3. actual global ToE closure,
4. actual strict-core selector closure,
5. actual global `QW-2191` discharge.
