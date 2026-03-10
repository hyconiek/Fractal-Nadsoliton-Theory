# P283 Current Actual Component 2 Explicit Three-Branch Provider Frontier Probe

Status: `P283_EXECUTED_CURRENT_ACTUAL_COMPONENT_2_EXPLICIT_THREE_BRANCH_PROVIDER_FRONTIER_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Probe whether the already explicit three component-2 provider branches still
contain one honest entering move on the current repo state.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| explicit fractal branch remains in scope | YES | `N294` keeps the branch in lane below support |
| explicit fractal branch enters component 2 | NO | the branch is map-rule blocked |
| explicit preobserver branch remains in scope | YES | `N295` keeps the branch in lane below support |
| explicit preobserver branch enters component 2 | NO | the branch is pair-index nonentering |
| explicit residual-datum branch remains in scope | YES | `N302` keeps the branch in lane below support |
| explicit residual-datum branch enters component 2 | NO | the branch is object-support blocked |
| an entering move remains inside the already explicit three-branch set | NO | no further honest positive lift is available there |
| next honest positive move requires new provider class or new blocker-cut | YES | the frontier is now sharp |

## Probe result

`P283` returns:

```text
explicit known three component-2 branches present: yes
explicit known three component-2 branches entering: no
```

## Consequence

The strongest honest current repo reading is:

1. all three already explicit branches remain scientifically live enough to be
   named,
2. none of the three branches currently enters component 2,
3. so the next honest positive move must come from outside that already
   explicit three-branch set.
