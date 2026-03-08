# P244 Current Actual Legacy-to-Strict Kernel Amplitude Nonabsorption Component Witness Probe

Status: `P244_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COMPONENT_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Test whether the current repo exports one first actual component-level
obstruction witness for the `T16` amplitude layer, while still remaining below
full amplitude obstruction.

## Input

`P244` reads:

1. `N65`,
2. `N70`,
3. `N83`,
4. `N263`,
5. `F154`.

## Probe question

Does the current repo export:

```text
A_abs_nonbridge_component_witness_v1 :
  (K_legacy_ont, K_strict_gate)
    -> amplitude_nonabsorption_component_obstruction_tag_v1
```

such that:

1. the legacy Weinberg amplitude role remains nontransferred,
2. the strict-side candidate object is present,
3. no explicit role-equivalence verdict is exported,
4. the full claim-specific legacy Weinberg frontier is already closed
   negatively,
5. the bridge/nonbridge frontier remains undecided?

## Expected outcome

If the witness is honest, the strongest expected current statement is:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COMPONENT_WITNESS_BELOW_FULL_AMPLITUDE_OBSTRUCTION_AFTER_P244
```

## Hard limits

Passing `P244` does not mean:

1. full amplitude obstruction is discharged,
2. strengthened nonbridge is discharged,
3. bridge is ruled out,
4. branch selection is justified,
5. selector closure follows.
