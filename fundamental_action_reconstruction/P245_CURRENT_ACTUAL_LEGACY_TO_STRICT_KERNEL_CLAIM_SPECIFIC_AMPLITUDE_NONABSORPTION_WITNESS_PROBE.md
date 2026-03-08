# P245 Current Actual Legacy-to-Strict Kernel Claim-Specific Amplitude Nonabsorption Witness Probe

Status: `P245_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_CLAIM_SPECIFIC_AMPLITUDE_NONABSORPTION_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

Test whether the current repo now exports one actual claim-specific amplitude
nonabsorption witness above the first component witness and still below full
amplitude obstruction.

## Input

`P245` reads:

1. `N264`,
2. `F155`.

## Probe question

Does the current repo export:

```text
A_abs_nonbridge_claim_specific_actual_witness_v1 :
  (K_legacy_ont, K_strict_gate)
    -> claim_specific_amplitude_nonabsorption_obstruction_target_v1
```

such that:

1. the first component witness is present,
2. the witness remains claim-specific,
3. the bridge/nonbridge frontier remains undecided,
4. full amplitude obstruction is still not discharged?

## Expected outcome

If the packet is honest, the strongest expected current statement is:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_CLAIM_SPECIFIC_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_WITNESS_BELOW_FULL_AMPLITUDE_OBSTRUCTION_AFTER_P245
```

## Hard limits

Passing `P245` does not mean:

1. full amplitude obstruction is discharged,
2. strengthened nonbridge is discharged,
3. bridge is ruled out,
4. branch selection is justified,
5. selector closure follows.
