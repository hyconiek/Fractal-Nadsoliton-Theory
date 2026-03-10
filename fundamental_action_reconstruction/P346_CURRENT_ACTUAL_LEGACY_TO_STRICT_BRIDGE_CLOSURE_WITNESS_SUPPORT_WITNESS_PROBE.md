# P346 Current Actual Legacy-To-Strict Bridge Closure Witness Support Witness Probe

Status: `P346_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_BRIDGE_CLOSURE_WITNESS_SUPPORT_WITNESS_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P346` checks whether the current repo exports one actual support witness
below the bridge-closure support packet while keeping the result:

1. below actual bridge discharge,
2. below branch selection,
3. below role transfer,
4. below strict-core selector closure,
5. below ToE closure.

## Fixed object under test

```text
Lambda_legacy_strict_bridge_closure_witness_support_witness_v1 :=
(
  Kappa_legacy_strict_bridge_closure_witness_support_packet_v1,
  B_legacy_strict_bridge_target_v1,
  Gamma_bridge_components_target_v1
)
```

## What P346 checks

`P346` checks only:

1. the support witness is explicitly exported,
2. it preserves the bridge/nonbridge bifurcation context,
3. it preserves split separation between `K_legacy_ont` and `K_strict_gate`,
4. it preserves the `N269` discipline that bridge is not a mandatory `T14`
   closure gate,
5. it remains below actual bridge derivation,
6. it remains below branch selection,
7. it remains below role transfer,
8. it remains below theory closure.

## Result

`P346` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_BRIDGE_CLOSURE_WITNESS_SUPPORT_WITNESS_BELOW_ACTUAL_BRIDGE_AFTER_P346
```

This means:

1. the positive bridge branch now has one actual support witness below the
   support packet and below the future-only closure target,
2. but it still does not export actual bridge discharge,
3. it still does not choose bridge over nonbridge,
4. it still does not justify strict-core selector closure or ToE closure.

## Hard limits

`P346` does not establish:

1. actual legacy-to-strict bridge derivation,
2. actual bridge-closure witness,
3. branch selection in favor of bridge,
4. legacy physical-role transfer,
5. strict-core selector closure,
6. global `QW-2191` discharge,
7. ToE closure.
