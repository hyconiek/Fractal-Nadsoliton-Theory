# P248 Current Actual Legacy-to-Strict Kernel Damping Nonrenormalization Obstruction Witness Probe

Status: `P248_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_DAMPING_NONRENORMALIZATION_OBSTRUCTION_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Test whether the current repo now exports one actual damping-layer obstruction
for the `T16` route.

## Input

`P248` reads:

1. `P47/N50`,
2. `N117`,
3. `N267`,
4. `F158`.

## Probe question

Does the current repo export:

```text
R_damp_nonbridge_actual_obstruction_witness_v1 :
  (K_legacy_ont, K_strict_gate)
    -> R_damp_nonbridge_obstruction_target_v1
```

such that:

1. no explicit `beta_tors -> (beta, eta)` translation rule is exported,
2. package-level nontransfer remains discharged,
3. amplitude obstruction is already discharged,
4. phase remains undischarged?

## Expected outcome

If the witness is honest, the strongest expected current statement is:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_DAMPING_NONRENORMALIZATION_OBSTRUCTION_WITNESS_AFTER_P248
```

## Hard limits

Passing `P248` does not mean:

1. phase obstruction is discharged,
2. strengthened nonbridge is discharged,
3. bridge is ruled out forever,
4. branch selection is justified,
5. selector closure follows.
