# P260 Current Actual Strict-Side Selector Ingredient Second Clause Extension Lift Probe

Status: `P260_EXECUTED_CURRENT_ACTUAL_STRICT_SIDE_SELECTOR_INGREDIENT_SECOND_CLAUSE_EXTENSION_LIFT_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

`P260` tests whether the current repo really exports the second-clause
extension-lift packet from `F169`, while keeping the result:

1. below actual `E_orient`,
2. below admissible `S_sel_int`,
3. below strict-core selector closure,
4. below ToE closure.

## What P260 checks

`P260` checks only:

1. the original strict-core second clause remains undischarged,
2. the strict-side admissibility principle is now accepted in
   `strict_extension_only` scope,
3. the first clause already has one explicit extension lift,
4. `S_sel_int_candidate_seed_v0` is now explicitly lifted only to one
   extension-scoped carrier-typed precursor status,
5. no claim identifies that lift with actual `E_orient`,
6. no claim identifies that lift with strict-core admission,
7. no strict-core selector closure is claimed.

## Result

`P260` returns:

```text
CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_SIDE_SELECTOR_INGREDIENT_SECOND_CLAUSE_EXTENSION_LIFT_PACKET_BELOW_E_ORIENT_AND_STRICT_ADMISSION_AFTER_P260
```

This means:

1. the second clause is no longer only one undischarged strict-core contract
   requirement,
2. it now also has one explicit extension-scoped carrier-typed precursor lift,
3. the repo still does not justify actual `E_orient`, strict-core admission,
   or closure language.

## Hard limits

`P260` does not establish:

1. actual `E_orient`,
2. admissible `S_sel_int`,
3. actual strict-core selector closure,
4. actual global `QW-2191` discharge,
5. actual ToE closure.
