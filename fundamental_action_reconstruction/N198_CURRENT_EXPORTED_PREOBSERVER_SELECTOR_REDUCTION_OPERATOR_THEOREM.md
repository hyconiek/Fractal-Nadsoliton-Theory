# N198 Current Exported Preobserver Selector Reduction Operator Theorem

Status: `N198_DISCHARGED_CURRENT_EXPORTED_PREOBSERVER_SELECTOR_REDUCTION_OPERATOR_THEOREM_NO_FALSE_PASS`
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

and its admissible selector bridge operator

```text
B_sel_preLM_v1
```

already export one admissible strict-core preobserver selector reduction
operator:

```text
R_sel_preLM_v1 : V_topo ⊕ L_int ⊕ M_int -> Q_sel_v1
```

## Scope

This theorem is scoped only to:

- current repo state,
- `S_preLM_strict_core_source_object_v1`,
- `E_orient_preLM_v1`,
- `B_sel_preLM_v1`,
- `R_sel_preLM_v1`,
- the strict preobserver carriers `[u_T, u_L, u_M]` and `[q_+, q_-]`.

## Proof basis

`P178` discharges the required reduction clauses:

1. the orientation input is admissible,
2. the selector bridge input is admissible,
3. the reduction map is derived only from those inputs,
4. it remains strict-core only,
5. it remains preobserver only,
6. it exports an internal selector sector,
7. it yields a positive plus-channel source response,
8. it yields a vanishing minus-channel source response,
9. it is bridge-ready for a later `O_sel`,
10. it is kernel-split-safe.

Therefore `R_sel_preLM_v1` is an admissible strict-core preobserver selector
reduction operator on current repo state.

## Hard limits

This theorem does not claim:

- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.
