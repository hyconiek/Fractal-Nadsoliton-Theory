# O3 A1_ext Coefficient Evaluation Rule

## Goal

Attach the first explicit coefficient-evaluation rule to the persisted operator
instance

- `A_1_ext(pair1)`

created in `O2`, without pretending that any coefficient values are already
numerically known.

## Starting Point

From `O2`:

- `A_1_ext(pair1) = [[a_1, b_1], [b_1, d_1]]`
- basis ordering is fixed as `[c_1, s_1]`
- the operator instance is persisted on the lane:
  - `hypothesis_extension_only`

From `H20`:

- coefficient semantics is already fixed.

What is still missing is the direct rule that says how the persisted matrix
instance is to be read when a future concrete export or evaluation becomes
available.

## Minimal Evaluation Rule

For the persisted `A_1_ext(pair1)` matrix written in the basis

- `[c_1, s_1]`

with entries

- `[[a_1, b_1], [b_1, d_1]]`

we define the first direct readout rule:

1. `eval(a_1) := read_entry(A_1_ext, 0, 0)`
2. `eval(b_1) := read_entry(A_1_ext, 0, 1)`
3. `eval(d_1) := read_entry(A_1_ext, 1, 1)`

Derived readouts are then:

4. `eval(tr(A_1)) := eval(a_1) + eval(d_1)`
5. `eval(Delta_1) := (eval(a_1) - eval(d_1), eval(b_1))`

## Scope Discipline

This rule does **not** assert that:

- the matrix entries are numeric,
- the matrix entries are exported from an operator chain,
- the selector-breaking test has succeeded,
- the operator belongs to strict core.

It only specifies how a future concrete `A_1_ext` instance is to be read.

## Result

The repository now contains the first explicit coefficient-evaluation rule for
`A_1_ext(pair1)`.

So the problem is no longer:

- missing readout rule for the persisted operator instance

but instead:

- the persisted entries remain symbolic placeholders, so the rule still cannot
  produce actual coefficient values.

## Honest Frontier

- `O3_B1 := the persisted A_1_ext instance now has an explicit coefficient-evaluation rule, but its entries remain symbolic placeholders, so no actual values for a_1, b_1, d_1, tr(A_1), or Delta_1 are produced yet`
- `O2_B1 := a persisted exported_composite_A_1 instance for A_1_ext on pair1 now exists, but its coefficient entries remain symbolic and unevaluated, so no selector-breaking test can yet be executed`
- `H28_B1 := the current repository state contains no computable operator-level source from which a_1, b_1, d_1 can be actually exported or evaluated for pair1, even though Route A provenance and coefficient semantics are already in place`
- `T12_B1 := the typing judgment with totality and uniqueness is specified but not discharged for the current selector track`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and equivalence/export map remain absent`
- `C32_B2 := raw cross-pair overlap route is formally degenerate under the strict orthonormal-disjoint mode scaffold and thus does not export alpha_12`

## Negative Claims Maintained

- no theorem-level PASS
- no full-closure PASS
- no claim that any coefficient value is numerically known
- no claim that the selector-breaking test has been executed
- no claim that `A_1_ext` belongs to strict core
