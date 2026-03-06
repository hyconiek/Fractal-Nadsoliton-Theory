# O4 A1_ext Entry Population Rule

## Goal

Specify the first explicit rule saying what counts as an admissible populated
entry for the persisted operator instance

- `A_1_ext(pair1) = [[a_1, b_1], [b_1, d_1]]`

without pretending that any such populated entry already exists.

## Starting Point

From `O2`:

- the persisted operator object `A_1_ext(pair1)` exists,
- its basis is fixed as `[c_1, s_1]`,
- its entries remain symbolic placeholders.

From `O3`:

- the readout rule for `a_1`, `b_1`, `d_1` is already explicit.

What is still missing is the criterion that tells us when one of these
placeholders becomes a legitimate populated entry.

## Minimal Admissible Population Rule

A placeholder entry `x in {a_1, b_1, d_1}` is considered **populated** only if
at least one of the following two admissible witness forms exists.

### Route P1: direct matrix-entry export

A persisted record exists with all of:

- `object_id = A_1_ext`
- `pair = pair1`
- `basis_order = [c_1, s_1]`
- `entry = (i,j)` matching the slot of `x`
- `value = v_x`
- `operator_origin = exported_composite_A_1`
- `lane = hypothesis_extension_only`
- `base_kernel_contains_obs = false`
- `strict_core_reinterpretation = false`
- `selector_smuggling = false`

### Route P2: basis-scalar export

A persisted scalar witness exists with all of:

- `x = <basis_i, A_1_ext basis_j>`
- `basis_i, basis_j in {c_1, s_1}`
- `value = v_x`
- the same provenance discipline as Route P1

## Slot Map

For the fixed basis `[c_1, s_1]`:

- `a_1` corresponds to `(0,0)` and `<c_1, A_1_ext c_1>`
- `b_1` corresponds to `(0,1)` and `<c_1, A_1_ext s_1>`
- `d_1` corresponds to `(1,1)` and `<s_1, A_1_ext s_1>`

## Result

The repository now contains the first explicit admissible entry-population rule
for `A_1_ext(pair1)`.

So the problem is no longer:

- missing criterion for what counts as a legal populated entry

but instead:

- no actual Route P1 or Route P2 witness has yet been exported for any of the
  three entries.

## Honest Frontier

- `O4_B1 := the persisted A_1_ext instance now has an explicit admissible entry-population rule, but no actual Route P1 or Route P2 witness exists yet for a_1, b_1, or d_1`
- `O3_B1 := the persisted A_1_ext instance now has an explicit coefficient-evaluation rule, but its entries remain symbolic placeholders, so no actual values for a_1, b_1, d_1, tr(A_1), or Delta_1 are produced yet`
- `O2_B1 := a persisted exported_composite_A_1 instance for A_1_ext on pair1 now exists, but its coefficient entries remain symbolic and unevaluated, so no selector-breaking test can yet be executed`
- `H28_B1 := the current repository state contains no computable operator-level source from which a_1, b_1, d_1 can be actually exported or evaluated for pair1, even though Route A provenance and coefficient semantics are already in place`
- `T12_B1 := the typing judgment with totality and uniqueness is specified but not discharged for the current selector track`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and equivalence/export map remain absent`
- `C32_B2 := raw cross-pair overlap route is formally degenerate under the strict orthonormal-disjoint mode scaffold and thus does not export alpha_12`

## Negative Claims Maintained

- no theorem-level PASS
- no full-closure PASS
- no claim that any entry has been populated
- no claim that any coefficient value is numerically known
- no claim that selector-breaking has been tested
- no claim that `A_1_ext` belongs to strict core
