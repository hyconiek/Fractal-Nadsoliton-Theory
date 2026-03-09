# T21 Strict-Side Selector Ingredient First Clause Extension Lift Spec

Status: `T21_PACKET_READY_STRICT_SIDE_SELECTOR_INGREDIENT_FIRST_CLAUSE_EXTENSION_LIFT_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N278`, the next honest strict-side question is:

```text
can the first strict-side admissibility clause now be lifted
into one weaker extension-scoped precursor clause
without pretending that the strict-core clause itself is discharged?
```

`T21` does not discharge the original strict-core clause.

It does something narrower:

1. writes one packet-ready theorem spec for a first-clause extension lift,
2. keeps the lift inside `strict_extension_only` scope,
3. keeps admissible `S_sel_int`, strict-core selector closure, and ToE
   closure out.

## Extension-lifted clause

### Informal statement

If:

1. the strict-side admissibility principle is accepted in
   `strict_extension_only` scope,
2. the first strict-side clause remains undischarged in strict core,
3. the repo exports one actual clause-support packet and one actual
   source-topology support family,

then the current seed candidate may be lifted only to:

```text
one extension-admissible selector-ingredient precursor
```

and not to admissible `S_sel_int`.

### Formal target

```text
T21_StrictSideSelectorIngredient_FirstClause_ExtensionLift

Assume:
  A1. the first strict-core clause genuinely_new_strict_core_source_object_required
      remains undischarged in strict core;
  A2. the strict-side admissibility principle is accepted at theory level in
      strict_extension_only scope;
  A3. one actual first-clause support packet is exported;
  A4. one actual source-topology selector support family is exported upstream
      of observer;
  A5. no theorem identifies S_sel_int_candidate_seed_v0 with admissible
      S_sel_int.

Then:
  C1. S_sel_int_candidate_seed_v0 may be lifted only to one extension-scoped
      selector-ingredient precursor clause target;
  C2. that lift remains weaker than admissible S_sel_int;
  C3. that lift remains below strict-core selector closure and ToE closure.
```

## Acceptance skeleton

This spec is acceptable only if all of the following stay explicit:

1. the original strict-core clause remains undischarged,
2. the lift is only extension-scoped,
3. no claim is made that `S_sel_int_candidate_seed_v0` is a genuine new
   strict-core source object,
4. no claim is made that `E_orient` already exists,
5. no claim is made that strict-core selector closure or ToE closure follows.

## Recommended next move

The correct next move is:

1. package one explicit extension-lift packet for the first clause,
2. test whether that packet is really exported,
3. keep the move below admissible `S_sel_int`.
