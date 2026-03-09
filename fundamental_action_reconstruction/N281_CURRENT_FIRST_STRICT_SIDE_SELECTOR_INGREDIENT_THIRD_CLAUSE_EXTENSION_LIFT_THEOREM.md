# N281 Current First Strict-Side Selector Ingredient Third Clause Extension Lift Theorem

Status: `N281_DISCHARGED_CURRENT_FIRST_STRICT_SIDE_SELECTOR_INGREDIENT_THIRD_CLAUSE_EXTENSION_LIFT_THEOREM_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Package theorem-level the strongest honest direct clause result available after
`P261`.

## Theorem-level conclusion

From `P261`, the current repo exports one actual extension-lift witness:

```text
Eta_strict_selector_clause3_extension_lift_actual_witness_v1 :
  S_sel_int_candidate_seed_v0
    -> strict_extension_selector_ingredient_precursor_clause3_target_v1
```

with the following scoped meaning:

1. the original strict-core third clause remains undischarged,
2. the accepted extension principle now permits one weaker source-seed-only
   precursor lift of that clause,
3. `S_sel_int_candidate_seed_v0` may therefore be treated as one
   extension-scoped source-seed-only precursor for later work,
4. the lift remains weaker than actual `E_orient`,
5. the lift remains weaker than actual `B_sel`, `R_sel`, and `O_sel`,
6. the lift remains weaker than admissible `S_sel_int`,
7. the lift remains below strict-core selector closure and ToE closure.

## What N281 does prove

`N281` proves only this narrower statement:

1. the repo now exports one actual direct extension-lift of the third
   strict-side admissibility clause,
2. this is stronger than leaving that clause only at contract level plus
   theory-level principle acceptance,
3. the new result still remains extension-only and non-closure.

## What N281 does not prove

`N281` does not prove:

1. the original strict-core third clause,
2. actual `E_orient`,
3. actual `B_sel`, `R_sel`, or `O_sel`,
4. admissible `S_sel_int`,
5. actual strict-core selector closure,
6. actual global selector closure,
7. actual global `QW-2191` discharge,
8. actual ToE closure,
9. legacy-to-strict bridge derivation.

## Consequence

The strongest honest current reading after `N281` is:

1. the first three minimal admissibility clauses now each have:
   - one undischarged strict-core requirement,
   - one extension-scoped precursor lift,
2. the repo still does not export the missing strict-core selector ingredient,
3. four original admissibility clauses remain unlifted,
4. among those four, the next original clause
   `strict-core only`
   is structurally nontrivial on the present `strict_extension_only` lane.
