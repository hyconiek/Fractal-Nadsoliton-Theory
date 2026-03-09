# N279 Current First Strict-Side Selector Ingredient First Clause Extension Lift Theorem

Status: `N279_DISCHARGED_CURRENT_FIRST_STRICT_SIDE_SELECTOR_INGREDIENT_FIRST_CLAUSE_EXTENSION_LIFT_THEOREM_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Package theorem-level the strongest honest direct clause result available after
`P259`.

## Theorem-level conclusion

From `P259`, the current repo exports one actual extension-lift witness:

```text
Psi_strict_selector_clause1_extension_lift_actual_witness_v1 :
  S_sel_int_candidate_seed_v0
    -> strict_extension_selector_ingredient_precursor_clause1_target_v1
```

with the following scoped meaning:

1. the original strict-core first clause remains undischarged,
2. the accepted extension principle now permits one weaker precursor lift of
   that clause,
3. `S_sel_int_candidate_seed_v0` may therefore be treated as one
   extension-scoped selector-ingredient precursor,
4. the lift remains weaker than admissible `S_sel_int`,
5. the lift remains below strict-core selector closure and ToE closure.

## What N279 does prove

`N279` proves only this narrower statement:

1. the repo now exports one actual direct extension-lift of the first
   strict-side admissibility clause,
2. this is stronger than leaving that clause only at obstruction level plus
   theory-level principle acceptance,
3. the new result still remains extension-only and non-closure.

## What N279 does not prove

`N279` does not prove:

1. the original strict-core first clause,
2. admissible `S_sel_int`,
3. actual strict-core selector closure,
4. actual global selector closure,
5. actual global `QW-2191` discharge,
6. actual ToE closure,
7. legacy-to-strict bridge derivation.

## Consequence

The strongest honest current reading after `N279` is:

1. the first clause now has both:
   - one strict-core obstruction,
   - and one extension-scoped precursor lift,
2. the repo still does not export the missing strict-core selector ingredient,
3. the next honest move is either a next-clause extension lift or a genuinely
   new strict-core ingredient attempt.
