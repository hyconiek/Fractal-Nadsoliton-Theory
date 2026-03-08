# P247 Current Actual Legacy-to-Strict Kernel Full Amplitude Nonabsorption Obstruction Witness Probe

Status: `P247_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_FULL_AMPLITUDE_NONABSORPTION_OBSTRUCTION_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Test whether the current repo now exports one actual full amplitude-layer
obstruction witness for the `T16` route.

## Input

`P247` reads:

1. `N50`,
2. `N117`,
3. `N266`,
4. `F157`.

## Probe question

Does the current repo export:

```text
A_abs_nonbridge_actual_obstruction_witness_v1 :
  (K_legacy_ont, K_strict_gate)
    -> A_abs_nonbridge_obstruction_target_v1
```

such that:

1. no exported explicit amplitude absorption map for `alpha_geo` is present,
2. package-level nontransfer remains discharged,
3. amplitude coverage over the current legacy role package is exported,
4. damping and phase layers remain undischarged?

## Expected outcome

If the witness is honest, the strongest expected current statement is:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_FULL_AMPLITUDE_NONABSORPTION_OBSTRUCTION_WITNESS_AFTER_P247
```

## Hard limits

Passing `P247` does not mean:

1. strengthened nonbridge is discharged,
2. bridge is ruled out forever,
3. branch selection is justified,
4. selector closure follows.
