# H16 Operator-Origin Partial Witness Audit

## Goal

Check whether either of the two admissible `operator_origin` values identified in `H13` already has any partial witness in the repository, under the `H15` separation between existing kernel feedback and the new `K_obs` extension lane.

## Admissible Values From H13

1. `exported_composite_A_1`
2. `pullback_from_E_1_G_light_R_mat_O_obs`

## Witness Standard

A partial witness may be any of the following:

- an explicit composite formula,
- a named carrier schema,
- a factorized route skeleton,
- a persistent candidate object,
- or a provenance record with at least some fields populated.

A partial witness is **not** enough for:

- provenance-valid `Route A`,
- selector discharge,
- strict-core promotion.

## Audit Result

### 1. `exported_composite_A_1`

This value has a **formula-level partial witness** in the repository.

Explicit occurrences already exist:

- `A_i = P_i E_i^* G_light^* R_mat^* O_obs R_mat G_light E_i P_i`
- specialized for pair1 as `A_1 = P_1 E_1^* G_light^* R_mat^* O_obs R_mat G_light E_1 P_1`
- minimal persisted candidate `A_1_cand = [[a_1,b_1],[b_1,d_1]]`
- partial provenance record for `A_1_cand`

Therefore `exported_composite_A_1` is **not empty**. It has a composite-formula and candidate-object witness.

However, it still lacks:

- an explicit exported matrix representative,
- evaluated coefficients `(a_1,b_1,d_1)`,
- provenance-valid operator origin.

### 2. `pullback_from_E_1_G_light_R_mat_O_obs`

This value has only a **factor-chain slot-level partial witness**.

Present in repo:

- named chain slots `E_1`, `G_light`, `R_mat`, `O_obs`,
- route-B construction schema,
- composite pullback formula referring to that chain.

Missing:

- exported action tables,
- matrix representatives,
- instantiated chain carriers,
- concrete pullback object on pair1.

Therefore this value is weaker than the composite candidate lane.

## Comparative Conclusion

The repository contains partial witnesses for **both** admissible values, but they are asymmetric:

- `exported_composite_A_1` has a stronger composite-formula / candidate-object witness,
- `pullback_from_E_1_G_light_R_mat_O_obs` has only a weaker factor-chain slot witness.

Neither one yields a provenance-valid `Route A` instance for pair1.

## Honest Frontier

- `H16_B1 := both admissible operator_origin values now have partial witnesses, but only at unequal strength (composite-formula/candidate-object vs factor-chain-slot), and neither witness reaches a provenance-valid Route A instance for pair1`
- `H15_B1 := existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback`
- `T12_B1 := the typing judgment with totality and uniqueness is specified but not discharged for the current selector track`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and sign-only export-map object exist, but target-slot population (theta_1,theta_2) remains absent`
- `C32_B2 := raw cross-pair overlap route is formally degenerate under the strict orthonormal-disjoint mode scaffold and thus does not export alpha_12`

## Negative Claims Maintained

- no theorem-level PASS
- no full-closure PASS
- no claim that either admissible `operator_origin` is instantiated
- no claim that the composite candidate is provenance-valid
- no claim that the factor-chain witness is computationally usable
- no claim that `QW-2191` is discharged
