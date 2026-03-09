# N280 Current First Strict-Side Selector Ingredient Second Clause Extension Lift Theorem

Status: `N280_DISCHARGED_CURRENT_FIRST_STRICT_SIDE_SELECTOR_INGREDIENT_SECOND_CLAUSE_EXTENSION_LIFT_THEOREM_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Package theorem-level the strongest honest direct clause result available after
`P260`.

## Theorem-level conclusion

From `P260`, the current repo exports one actual extension-lift witness:

```text
Chi_strict_selector_clause2_extension_lift_actual_witness_v1 :
  S_sel_int_candidate_seed_v0
    -> strict_extension_selector_ingredient_precursor_clause2_target_v1
```

with the following scoped meaning:

1. the original strict-core second clause remains undischarged,
2. the accepted extension principle now permits one weaker carrier-typed
   precursor lift of that clause,
3. `S_sel_int_candidate_seed_v0` may therefore be treated as one
   extension-scoped carrier-typed precursor for later export work,
4. the lift remains weaker than actual `E_orient`,
5. the lift remains weaker than admissible `S_sel_int`,
6. the lift remains below strict-core selector closure and ToE closure.

## What N280 does prove

`N280` proves only this narrower statement:

1. the repo now exports one actual direct extension-lift of the second
   strict-side admissibility clause,
2. this is stronger than leaving that clause only at contract level plus
   theory-level principle acceptance,
3. the new result still remains extension-only and non-closure.

## What N280 does not prove

`N280` does not prove:

1. the original strict-core second clause,
2. actual `E_orient`,
3. admissible `S_sel_int`,
4. actual strict-core selector closure,
5. actual global selector closure,
6. actual global `QW-2191` discharge,
7. actual ToE closure,
8. legacy-to-strict bridge derivation.

## Consequence

The strongest honest current reading after `N280` is:

1. the first two minimal admissibility clauses now each have:
   - one undischarged strict-core requirement,
   - one extension-scoped precursor lift,
2. the repo still does not export the missing strict-core selector ingredient,
3. the next honest move is either a next-clause extension lift or a genuinely
   new strict-core ingredient attempt.
