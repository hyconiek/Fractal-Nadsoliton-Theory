# N200 Current Preobserver To Emergent Observer Coarse Graining Theorem

Status: `N200_DISCHARGED_CURRENT_PREOBSERVER_TO_EMERGENT_OBSERVER_COARSE_GRAINING_THEOREM_NO_FALSE_PASS`
As of: `2026-03-08`

## Theorem statement

Within the current repo state,

```text
S_preLM_strict_core_source_object_v1
```

its admissible orientation datum

```text
E_orient_preLM_v1
```

its admissible selector bridge operator

```text
B_sel_preLM_v1
```

its admissible selector reduction operator

```text
R_sel_preLM_v1
```

its admissible selector output operator

```text
O_sel_preLM_v1
```

already export one admissible strict-core preobserver-to-emergent-observer
coarse-graining operator:

```text
C_obs_limit_preLM_v1 : Q_out_v1 -> Y_obs_limit_v1
```

## Scope

This theorem is scoped only to:

- current repo state,
- `S_preLM_strict_core_source_object_v1`,
- `E_orient_preLM_v1`,
- `B_sel_preLM_v1`,
- `R_sel_preLM_v1`,
- `O_sel_preLM_v1`,
- `C_obs_limit_preLM_v1`,
- the strict preobserver carriers `[u_T, u_L, u_M]`, `[q_+, q_-]`,
  `[o_+, o_-]`, and the observer-limit carrier `[y_bias, y_total]`.

## Proof basis

`P180` discharges the required coarse-graining clauses:

1. the selector-output input is admissible,
2. the coarse-graining map is derived only from that output,
3. it remains strict-core only,
4. it remains preobserver-to-observer-limit only,
5. it exports an observer-limit coarse-grained sector,
6. it yields a positive macroscopic bias readout,
7. it yields a positive total-signal readout,
8. it keeps observer information deficit downstream,
9. it is kernel-split-safe.

Therefore `C_obs_limit_preLM_v1` is an admissible strict-core preobserver to
emergent-observer coarse-graining operator on current repo state.

## Hard limits

This theorem does not claim:

- actual emergent observer construction,
- selector closure,
- `QW-2191` discharge,
- ToE closure.
