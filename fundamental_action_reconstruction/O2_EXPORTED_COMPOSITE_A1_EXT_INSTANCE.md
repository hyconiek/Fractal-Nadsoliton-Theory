# O2 Exported Composite A1_ext Instance

## Goal

Create the first persisted instance of the missing operator object

- `A_1_ext`

in the admissible mode:

- `exported_composite_A_1`

without pretending that coefficient values are already known.

## Scope

This step does **not** compute:

- `a_1`
- `b_1`
- `d_1`

It only materializes the operator instance in the currently strongest admissible
mode, using the already established provenance and coefficient semantics.

## Carrier

- pair: `pair1 = (c_1, s_1)`
- plane: `V_1 = span{c_1, s_1}`
- basis ordering: `[c_1, s_1]`

## Admissible Mode Chosen

- `operator_origin = exported_composite_A_1`

This choice is justified because:

1. it is one of the two admissible origins from `H13`,
2. it has the stronger witness from `H16/H17`,
3. it already reached provenance-valid `Route A` status in `H18`.

## Persisted Instance

The persisted operator instance is:

- `A_1_ext(pair1) = [[a_1, b_1], [b_1, d_1]]`

with:

- coefficient semantics inherited from `H20`,
- provenance inherited from `H18`,
- lane discipline:
  - `hypothesis_extension_only`

## Result

The repository now contains a persisted operator instance for `A_1_ext` in the
`exported_composite_A_1` mode.

This is the first step that turns the missing operator from a pure spec into an
actual object.

## Honest Frontier

- `O2_B1 := a persisted exported_composite_A_1 instance for A_1_ext on pair1 now exists, but its coefficient entries remain symbolic and unevaluated, so no selector-breaking test can yet be executed`
- `H28_B1 := the current repository state contains no computable operator-level source from which a_1, b_1, d_1 can be actually exported or evaluated for pair1, even though Route A provenance and coefficient semantics are already in place`
- `T12_B1 := the typing judgment with totality and uniqueness is specified but not discharged for the current selector track`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and equivalence/export map remain absent`
- `C32_B2 := raw cross-pair overlap route is formally degenerate under the strict orthonormal-disjoint mode scaffold and thus does not export alpha_12`

## Negative Claims Maintained

- no theorem-level PASS
- no full-closure PASS
- no claim that `a_1`, `b_1`, `d_1` have values
- no claim that `A_1_ext` is strict core
- no claim that any anisotropy test has been completed
