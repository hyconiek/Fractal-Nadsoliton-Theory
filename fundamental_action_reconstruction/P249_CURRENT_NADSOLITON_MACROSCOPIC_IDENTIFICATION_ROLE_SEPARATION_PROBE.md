# P249 Current Nadsoliton Macroscopic Identification Role Separation Probe

Status: `P249_EXECUTED_CURRENT_NADSOLITON_MACROSCOPIC_IDENTIFICATION_ROLE_SEPARATION_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Test whether the current repo now exports one actual role-separation packet
that reclassifies the `T15/T16` deadlock without claiming kernel identity or
selector closure.

## Input

`P249` reads:

1. `N117`,
2. `N260`,
3. `N263`,
4. `N268`,
5. `T17`,
6. `F159`.

## Probe question

Does the current repo export:

```text
Omega_nad_role_separation_actual_witness_v1 :
  (K_legacy_ont, K_strict_gate)
    -> nadsoliton_macro_identification_role_separation_target_v1
```

such that:

1. `K_legacy_ont` is fixed only as a macroscopic identification tool,
2. `K_strict_gate` is fixed only as a strict source-topology kernel,
3. lack of cross-kernel absorption is not treated as a ToE-failure witness by
   itself,
4. bridge/nonbridge remains open but no longer acts as a mandatory `T14`
   closure gate,
5. no kernel identity or legacy physical-role transfer is claimed?

## Expected outcome

If the packet is honest, the strongest expected current statement is:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_NADSOLITON_MACROSCOPIC_IDENTIFICATION_ROLE_SEPARATION_PACKET_WITHDRAWING_THE_T15_T16_DEADLOCK_AS_A_MANDATORY_T14_CLOSURE_GATE_AFTER_P249
```

## Hard limits

Passing `P249` does not mean:

1. bridge is proved,
2. strengthened nonbridge is proved,
3. strict-core selector closure is proved,
4. global `QW-2191` is discharged,
5. ToE is closed.
