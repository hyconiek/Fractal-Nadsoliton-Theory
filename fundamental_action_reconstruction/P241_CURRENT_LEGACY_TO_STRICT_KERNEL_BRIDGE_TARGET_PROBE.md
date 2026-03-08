# P241 Current Legacy-To-Strict Kernel Bridge Target Probe

Status: `P241_EXECUTED_CURRENT_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

`P241` tests whether the current repo really exports one explicit future-only Legacy-To-Strict Kernel Bridge target, while keeping the result:

1. below actual bridge derivation,
2. below any legacy physical-role transfer theorem,
2. below global `QW-2191` discharge,
3. below strict-core selector closure.

## Fixed input

Input component packet:

```text
Gamma_bridge_components_target_v1 :=
(
  A_abs_margin_target_v1,
  R_damp_renorm_target_v1,
  P_conformal_shift_target_v1
)
```

Target relation under test:

```text
B_legacy_strict_bridge_target_v1 : K_legacy_ont -> K_strict_gate
```

## What P241 checks

`P241` checks only:

1. the packet is explicitly exported as a target,
2. the separation between `K_legacy_ont` and `K_strict_gate` is maintained without silent identification,
3. the packet remains future-only,
4. the non-bridge branch remains explicit,
5. no legacy physical-role transfer is claimed,
6. observer remains downstream of the bridge target,
7. the result remains below actual legacy-to-strict bridge derivation,
8. the result remains below global `QW-2191` discharge.

## Result

`P241` returns:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_BELOW_ACTUAL_BRIDGE_DISCHARGE_AFTER_P241
```

This means:

1. the current repo exports one future-only positive bridge-branch target,
2. but it still does not export an actual bridge derivation,
3. it still does not transfer legacy physical roles,
4. and it still does not justify selector closure or global `QW-2191`
   discharge.

## Hard limits

`P241` does not establish:

1. actual legacy-to-strict bridge derivation,
2. legacy physical-role transfer,
3. actual strict-core selector closure,
4. current global `QW-2191` discharge,
5. ToE closure.
