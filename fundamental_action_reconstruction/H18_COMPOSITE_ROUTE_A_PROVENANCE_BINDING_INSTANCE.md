# H18 Composite Route A Provenance Binding Instance

## Goal

Populate the single remaining provenance-binding step identified in `H17`, namely:

- bind `A_1_cand` to
- `operator_origin = exported_composite_A_1`

inside one provenance-valid `Route A` record for `pair1`, while keeping the lane explicitly outside strict core.

## Preconditions From Earlier Steps

Already established before this step:

1. `lane = hypothesis_extension_only`
2. `base_kernel_contains_obs = false`
3. `selector_smuggling = false`
4. `strict_core_reinterpretation = false`
5. `A_1_cand` exists as the current composite candidate on `pair1`
6. the stronger admissible `operator_origin` lane is `exported_composite_A_1`
7. no evaluated coefficient claim is made for `(a_1,b_1,d_1)`
8. no selector discharge claim is made

## Binding Performed

The provenance-valid Route A binding for `pair1` is now instantiated as:

- `operator_origin = exported_composite_A_1`
- `current_composite_export_witness = A_1_cand`
- `carrier_basis = {c_1, s_1}`
- `lane = hypothesis_extension_only`
- `strict_core_reinterpretation = false`
- `base_kernel_contains_obs = false`

## Result

A first provenance-valid `Route A` witness now exists for `pair1`, but only on the hypothesis-extension lane.

This witness does **not** imply:

- evaluated coefficients,
- nontrivial anisotropy,
- selector breaking,
- strict-core closure,
- theorem-level PASS.

It only means that the stronger composite witness is no longer blocked at the provenance-binding level.

## Honest Frontier

- `H18_B1 := a provenance-valid Route A witness now exists on the hypothesis-extension lane for pair1, but no evaluated coefficient triple (a_1,b_1,d_1) has yet been extracted from it, so no O(2)-breaking test has been executed`
- `H15_B1 := existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback`
- `T12_B1 := the typing judgment with totality and uniqueness is specified but not discharged for the current selector track`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and sign-only export-map object exist, but target-slot population (theta_1,theta_2) remains absent`
- `C32_B2 := raw cross-pair overlap route is formally degenerate under the strict orthonormal-disjoint mode scaffold and thus does not export alpha_12`

## Negative Claims Maintained

- no theorem-level PASS
- no full-closure PASS
- no claim that `A_1_cand` coefficients are evaluated
- no claim that selector degeneracy is discharged
- no claim that this witness belongs to strict core
- no claim that `QW-2191` is discharged
