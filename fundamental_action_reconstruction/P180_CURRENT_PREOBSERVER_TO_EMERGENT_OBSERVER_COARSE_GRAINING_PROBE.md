# P180 Current Preobserver To Emergent Observer Coarse Graining Probe

Status: `CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_PREOBSERVER_TO_EMERGENT_OBSERVER_COARSE_GRAINING_OPERATOR_FROM_O_SEL_PRELM_V1_AFTER_P180`
As of: `2026-03-08`

## Goal

Test whether the current repo now exports one actual admissible strict-core
preobserver-to-emergent-observer coarse-graining operator:

```text
O_sel_preLM_v1 -> C_obs_limit_preLM_v1
```

without pretending that an actual emergent observer, selector closure, or
`QW-2191` discharge already exist.

## Result

The probe passes on current repo state.

The repo now exports:

```text
C_obs_limit_preLM_v1 : Q_out_v1 -> Y_obs_limit_v1
```

with ordered observer-limit basis:

```text
y_bias, y_total
```

For the current source object:

```text
C_obs_limit_preLM_v1(
  O_sel_preLM_v1(R_sel_preLM_v1(S_preLM_strict_core_source_object_v1))
)
  = y_bias_v1 y_bias + y_total_v1 y_total
```

where:

```text
y_bias_v1  > 0
y_total_v1 > 0
```

## What this proves

The current repo exports one actual downstream coarse-graining map that:

1. is derived only from the admissible selector-output operator,
2. remains strict-core only,
3. remains preobserver-to-observer-limit only,
4. exports an observer-limit coarse-grained sector,
5. gives a positive macroscopic bias readout,
6. gives a positive total-signal readout,
7. keeps observer information deficit downstream,
8. is kernel-split-safe.

## Hard limits

This probe does not claim:

- actual emergent observer construction,
- selector closure,
- `QW-2191` discharge,
- ToE closure.
