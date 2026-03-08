# P243 Current Actual Legacy-to-Strict Kernel Bifurcated Frontier Probe

Status: `P243_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Test whether the current repo now exports one actual bifurcated frontier
packet for the highest-priority `legacy -> strict bridge or non-bridge`
question, while still refusing a false branch winner.

## Input

`P243` reads:

1. `N261` for the future-only positive bridge branch,
2. `N262` for the future-only negative nonbridge-strengthening branch,
3. `F153` for the actual frontier packet that bundles them.

## Probe question

Does the current repo export:

```text
Xi_legacy_strict_frontier_bifurcation_packet_v1
```

such that:

1. the bridge branch is explicit,
2. the nonbridge-strengthening branch is explicit,
3. both branches remain future-only,
4. neither branch is currently discharged,
5. no current branch-selection theorem is exported?

## Expected outcome

If the packet is honest, the strongest expected current statement is:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PACKET_WITH_NO_JUSTIFIED_BRANCH_SELECTION_AFTER_P243
```

## Hard limits

Passing `P243` does not mean:

1. bridge is proved,
2. strengthened nonbridge is proved,
3. one branch now wins,
4. selector closure follows,
5. global `QW-2191` is discharged.
