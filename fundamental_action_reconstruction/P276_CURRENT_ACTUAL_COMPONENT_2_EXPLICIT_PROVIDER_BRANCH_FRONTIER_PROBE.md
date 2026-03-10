# P276 Current Actual Component 2 Explicit Provider Branch Frontier Probe

Status: `P276_EXECUTED_CURRENT_ACTUAL_COMPONENT_2_EXPLICIT_PROVIDER_BRANCH_FRONTIER_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Probe whether the already explicit component-2 provider branches still contain
one honest entering move on the current repo state.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| explicit fractal branch remains in scope | YES | `N294` keeps the branch in lane below support |
| explicit fractal branch enters component 2 | NO | the branch is map-rule blocked |
| explicit preobserver branch remains in scope | YES | `N295` keeps the branch in lane below support |
| explicit preobserver branch enters component 2 | NO | the branch is pair-index nonentering |
| an entering move remains inside the already explicit two-branch set | NO | no further honest positive lift is available there |
| next honest positive move requires third provider class or new blocker-cut | YES | the frontier is now sharp |

## Probe result

`P276` returns:

```text
explicit known component-2 branches present: yes
explicit known component-2 branches entering: no
```

## Consequence

The strongest honest current repo reading is:

1. both already explicit branches remain scientifically live enough to be
   named,
2. neither branch currently enters component 2,
3. so the next honest positive move must come from outside that already
   explicit two-branch set.
