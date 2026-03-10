# T99 Current Legacy-To-Strict Bridge Closure Witness Target Specification

Status: `T99_CURRENT_LEGACY_TO_STRICT_BRIDGE_CLOSURE_WITNESS_TARGET_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N261`, the positive bridge branch is explicit only as a future-only
kernel-comparison target.

After `N263`, the bridge/nonbridge frontier is explicit on both sides.

After `N269`, bridge/nonbridge is no longer treated as a mandatory closure
gate for the `T14` selector lane.

Therefore the next honest bridge-side move is not:

1. actual bridge discharge,
2. selector closure,
3. global `QW-2191` discharge,
4. ToE closure.

It is only:

```text
freeze one explicit future-only bridge-closure witness target
for the positive bridge branch itself
```

where "closure" means:

```text
closure of the positive bridge branch as a comparison-frontier route,
not closure of strict-core selector theory
and not closure of ToE
```

## Fixed target

Reuse:

```text
B_legacy_strict_bridge_target_v1 : K_legacy_ont -> K_strict_gate
```

and:

```text
Gamma_bridge_components_target_v1 :=
(
  A_abs_margin_target_v1,
  R_damp_renorm_target_v1,
  P_conformal_shift_target_v1
)
```

Freeze one explicit future-only target:

```text
Omega_legacy_strict_bridge_closure_witness_target_v1 :
(
  B_legacy_strict_bridge_target_v1,
  Gamma_bridge_components_target_v1
)
  -> declared_scope_legacy_strict_bridge_closure_target_v1
```

## Meaning

`declared_scope_legacy_strict_bridge_closure_target_v1` means only:

1. the positive bridge branch would be resolved at kernel-comparison level,
2. the three target components would be jointly packaged into one declared
   bridge-resolution object,
3. the branch would remain explicitly comparison-frontier-scoped.

It does not mean:

1. legacy physical-role transfer,
2. strict-core selector closure,
3. global selector closure,
4. global `QW-2191` discharge,
5. ToE closure.

## Scope discipline

`T99` remains guardrail-safe only if all of the following stay explicit:

1. `K_legacy_ont` and `K_strict_gate` remain split without silent
   identification,
2. the nonbridge branch remains open,
3. this branch is not promoted into a mandatory `T14` closure gate,
4. observer remains downstream,
5. any future role-transfer theorem remains separate.

## Hard limits

`T99` does not specify:

1. actual bridge derivation,
2. actual bridge-closure witness,
3. legacy physical-role transfer,
4. strict-core selector closure,
5. global `QW-2191` discharge,
6. ToE closure.
