# P344 Current Legacy-To-Strict Bridge Closure Witness Target Probe

Status: `P344_EXECUTED_CURRENT_LEGACY_TO_STRICT_BRIDGE_CLOSURE_WITNESS_TARGET_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P344` checks whether the current repo exports one explicit future-only
bridge-closure witness target while keeping the result:

1. below actual bridge discharge,
2. below branch selection,
3. below legacy physical-role transfer,
4. below strict-core selector closure,
5. below ToE closure.

## Fixed input

Input relation:

```text
B_legacy_strict_bridge_target_v1 : K_legacy_ont -> K_strict_gate
```

Input component packet:

```text
Gamma_bridge_components_target_v1 :=
(
  A_abs_margin_target_v1,
  R_damp_renorm_target_v1,
  P_conformal_shift_target_v1
)
```

Target under test:

```text
Omega_legacy_strict_bridge_closure_witness_target_v1 :
(
  B_legacy_strict_bridge_target_v1,
  Gamma_bridge_components_target_v1
)
  -> declared_scope_legacy_strict_bridge_closure_target_v1
```

## What P344 checks

`P344` checks only:

1. the new object is exported as a future-only target,
2. the bridge/nonbridge bifurcation remains explicit,
3. no branch-selection theorem is implied,
4. no silent kernel identification is implied,
5. no legacy physical-role transfer is claimed,
6. no strict-core selector closure is claimed,
7. `N269` remains explicit: the bridge route is not reintroduced as a
   mandatory `T14` closure gate,
8. the result remains comparison-frontier-scoped.

## Result

`P344` returns:

```text
CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_LEGACY_TO_STRICT_BRIDGE_CLOSURE_WITNESS_TARGET_BELOW_ACTUAL_BRIDGE_AND_BELOW_THEORY_CLOSURE_AFTER_P344
```

This means:

1. the positive bridge branch now has one explicit closure-target object,
2. but it still does not export actual bridge discharge,
3. it still does not select bridge over nonbridge,
4. it still does not justify selector closure or ToE closure.

## Hard limits

`P344` does not establish:

1. actual legacy-to-strict bridge derivation,
2. actual bridge-closure witness,
3. branch selection in favor of bridge,
4. legacy physical-role transfer,
5. strict-core selector closure,
6. global `QW-2191` discharge,
7. ToE closure.
