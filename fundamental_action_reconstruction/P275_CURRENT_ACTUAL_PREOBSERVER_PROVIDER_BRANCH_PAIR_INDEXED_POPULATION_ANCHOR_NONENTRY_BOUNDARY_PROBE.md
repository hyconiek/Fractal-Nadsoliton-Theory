# P275 Current Actual Preobserver Provider Branch Pair-Indexed Population Anchor Nonentry Boundary Probe

Status: `P275_EXECUTED_CURRENT_ACTUAL_PREOBSERVER_PROVIDER_BRANCH_PAIR_INDEXED_POPULATION_ANCHOR_NONENTRY_BOUNDARY_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Probe whether the existing preobserver provider branch currently enters the
pair-indexed population-anchor layer required by `T26` component 2.

## Probe table

| Check | Verdict | Meaning |
|---|---|---|
| preobserver provider branch exported | YES | `F73/F74/F75/F76` and `N179/N180/N181/N182` are in scope |
| branch stays upstream and observer-excluded | YES | current branch remains guardrail-consistent |
| branch exports pair-indexed reduction to `[pair1,pair2]` | NO | no pair-indexed reduction is exported on that branch |
| branch exports `theta_1/theta_2`-typed population layer | NO | no such reduction is exported on that branch |
| branch links itself to `R1/C48/C49` population scaffold | NO | no codomain-side population-anchor coupling is exported |
| branch exports actual component-2 support | NO | the branch remains below `pair_indexed_population_anchor_target_v1` |

## Probe result

`P275` returns:

```text
preobserver branch present: yes
preobserver branch enters component 2: no
```

## Consequence

The strongest honest current repo reading is:

1. the preobserver branch remains one valid future provider route,
2. but it is still nonentering for the pair-indexed population-anchor layer,
3. so it cannot yet be cited as actual support for component 2.
