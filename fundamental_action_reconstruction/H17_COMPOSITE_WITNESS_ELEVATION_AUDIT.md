# H17 Composite Witness Elevation Audit

## Goal

Determine whether the stronger partial witness for `operator_origin = exported_composite_A_1` can already be elevated to a provenance-valid `Route A` witness for `pair1`, without false reinterpretation of strict core.

## Starting Point

From `H16`:

- `exported_composite_A_1` has the stronger witness,
- because the repository already contains:
  - the composite formula for `A_1`,
  - the candidate carrier `A_1_cand`,
  - a partial provenance record,
- while the factor-chain lane is weaker and remains only slot-level.

## Elevation Test

A provenance-valid elevation from partial witness to `Route A` witness requires all of the following to be fixed simultaneously:

1. `lane = hypothesis_extension_only`
2. `base_kernel_contains_obs = false`
3. `operator_origin = exported_composite_A_1`
4. `selector_smuggling = false`
5. `strict_core_reinterpretation = false`
6. a direct statement that `A_1_cand` is the current composite export witness for pair1,
7. no claim that coefficients `(a_1,b_1,d_1)` are already evaluated,
8. no claim that this witness discharges selector degeneracy.

## Audit Result

The repository already satisfies conditions (1), (2), (4), (5), (7), and (8) at the level of packet-ready statements.

It does **not** yet satisfy condition (3) and (6) together in one provenance-valid record.

More precisely:

- `operator_origin` has been reduced to a finite admissible set,
- but it has not yet been populated as `exported_composite_A_1` inside a provenance-valid record,
- and `A_1_cand` has not yet been explicitly elevated from candidate object to current composite export witness for pair1 under the extension-only lane.

## Conclusion

The stronger composite witness can be elevated further than the factor-chain witness, but it is still blocked by one remaining binding step:

- explicit provenance-valid binding of `A_1_cand` to `operator_origin = exported_composite_A_1`

This is now the dominant unresolved issue on the `Route A` side.

## Honest Frontier

- `H17_B1 := the stronger composite witness for exported_composite_A_1 is one explicit provenance-binding step away from a provenance-valid Route A witness, but that binding has not yet been populated in the current record`
- `H15_B1 := existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback`
- `T12_B1 := the typing judgment with totality and uniqueness is specified but not discharged for the current selector track`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and sign-only export-map object exist, but target-slot population (theta_1,theta_2) remains absent`
- `C32_B2 := raw cross-pair overlap route is formally degenerate under the strict orthonormal-disjoint mode scaffold and thus does not export alpha_12`

## Negative Claims Maintained

- no theorem-level PASS
- no full-closure PASS
- no claim that `A_1_cand` is already provenance-valid
- no claim that the composite witness is strict-core native
- no claim that selector degeneracy is discharged
- no claim that `QW-2191` is discharged
