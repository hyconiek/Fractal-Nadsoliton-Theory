# N541 Current Second Admissibility Clause Discharge Theorem (S_sel_int strict‑core source object v1)

Status: `N541_DISCHARGED_CURRENT_SECOND_ADMISSIBILITY_CLAUSE_THEOREM_FOR_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_NO_FALSE_PASS`  
As of: `2026-03-16`

## Theorem statement

Within the current repo state,
`S_sel_int_strict_core_source_object_v1` discharges the second admissibility
clause:

```text
carrier_typed_enough_for_later_E_orient_export_required
```

and only that second clause (in addition to the already discharged first
clause).

## Scope

This theorem is scoped only to:

- current repo state,
- `S_sel_int_strict_core_source_object_v1`,
- the second admissibility clause of the `F34` minimal admissibility contract.

## Proof basis

1. `N540` discharges the first admissibility clause for `S_sel_int` by the
   existence of a genuinely new strict‑core exported source object.
2. `F649` freezes a minimal explicit typed carrier interface and a future
   `E_orient` target frame for that exported object, without exporting
   `E_orient`.
3. `P649` verifies that the second clause is satisfied in the narrow typed
   carrier sense: explicit carrier typing, explicit future target frame,
   nonzero sine support marker, and `E_orient` not falsely claimed.

Therefore the second clause is discharged on current repo state for this
exported object.

## Hard limits

This theorem does not claim:

- full admissibility of `S_sel_int`,
- actual `E_orient`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

