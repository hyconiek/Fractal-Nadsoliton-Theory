# P246 Current Actual Legacy-to-Strict Kernel Amplitude Nonabsorption Coverage Packet Probe

Status: `P246_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COVERAGE_PACKET_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Test whether the current repo now exports one actual amplitude-coverage packet
over the currently closed `alpha_geo`-bearing legacy role package, while still
remaining below full `A_abs` obstruction.

## Input

`P246` reads:

1. `N83`,
2. `N99`,
3. `N115`,
4. `N116`,
5. `N265`,
6. `F156`.

## Probe question

Does the current repo export:

```text
Kappa_abs_nonbridge_coverage_packet_v1
```

such that:

1. all three claim-specific legacy physical-role frontiers are closed
   negatively,
2. the legacy physical-role package is closed negatively,
3. one actual claim-specific amplitude witness is present,
4. full `A_abs` obstruction is still not discharged?

## Expected outcome

If the packet is honest, the strongest expected current statement is:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COVERAGE_PACKET_BELOW_FULL_A_ABS_OBSTRUCTION_AFTER_P246
```

## Hard limits

Passing `P246` does not mean:

1. full amplitude obstruction is discharged,
2. strengthened nonbridge is discharged,
3. bridge is ruled out,
4. branch selection is justified,
5. selector closure follows.
