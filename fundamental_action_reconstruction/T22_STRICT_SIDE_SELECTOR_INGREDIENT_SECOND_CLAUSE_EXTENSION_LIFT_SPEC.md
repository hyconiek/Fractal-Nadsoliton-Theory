# T22 Strict-Side Selector Ingredient Second Clause Extension Lift Spec

Status: `T22_PACKET_READY_STRICT_SIDE_SELECTOR_INGREDIENT_SECOND_CLAUSE_EXTENSION_LIFT_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N279`, the next honest strict-side question is:

```text
can the second strict-side admissibility clause now be lifted
into one weaker extension-scoped carrier-typed precursor clause
without pretending that E_orient or admissible S_sel_int already exists?
```

`T22` does not discharge the original strict-core second clause.

It does something narrower:

1. writes one packet-ready theorem spec for a second-clause extension lift,
2. keeps the lift inside `strict_extension_only` scope,
3. keeps `E_orient`, admissible `S_sel_int`, strict-core selector closure,
   and ToE closure out.

## Extension-lifted clause

### Informal statement

If:

1. the strict-side admissibility principle is accepted in
   `strict_extension_only` scope,
2. the second strict-side clause
   `carrier_typed_enough_for_later_export`
   remains undischarged in strict core,
3. the repo already exports one first-clause extension lift,
4. the frozen seed candidate still exists only as one minimal carrier pair,
5. one actual source-topology support family remains exported upstream of
   observer,

then the current seed candidate may be lifted only to:

```text
one extension-scoped carrier-typed precursor for later export
```

and not to actual `E_orient`, not to admissible `S_sel_int`.

### Formal target

```text
T22_StrictSideSelectorIngredient_SecondClause_ExtensionLift

Assume:
  A1. the second strict-core clause carrier_typed_enough_for_later_export
      remains undischarged in strict core;
  A2. the strict-side admissibility principle is accepted at theory level in
      strict_extension_only scope;
  A3. one actual first-clause extension lift is already exported;
  A4. S_sel_int_candidate_seed_v0 remains frozen only as the minimal carrier
      pair from F36;
  A5. one actual source-topology selector support family is exported upstream
      of observer;
  A6. no theorem exports E_orient from S_sel_int_candidate_seed_v0 and no
      theorem identifies that seed with admissible S_sel_int.

Then:
  C1. S_sel_int_candidate_seed_v0 may be lifted only to one extension-scoped
      carrier-typed precursor clause target for later export;
  C2. that lift remains weaker than actual E_orient export and weaker than
      admissible S_sel_int;
  C3. that lift remains below strict-core selector closure and ToE closure.
```

## Acceptance skeleton

This spec is acceptable only if all of the following stay explicit:

1. the original strict-core second clause remains undischarged,
2. the lift is only extension-scoped,
3. no claim is made that actual `E_orient` already exists,
4. no claim is made that downstream `B_sel`, `R_sel`, or `O_sel` already
   follow,
5. no claim is made that `S_sel_int_candidate_seed_v0` is already admissible
   `S_sel_int`,
6. no claim is made that strict-core selector closure or ToE closure follows.

## Recommended next move

The correct next move is:

1. package one explicit extension-lift packet for the second clause,
2. test whether that packet is really exported,
3. keep the move below actual `E_orient` and below admissible `S_sel_int`.
