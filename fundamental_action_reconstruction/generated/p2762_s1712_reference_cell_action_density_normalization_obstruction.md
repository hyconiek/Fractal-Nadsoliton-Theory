# P2762/S1712 reference-cell and action-density normalization obstruction

Status: `P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE`

## Finite scale-orbit witness
- tested_positive_reference_lengths=[0.25, 0.5, 1.0, 2.0, 4.0]
- distinct_dimensionalizations=5
- lambda_invariant_under_ell=True
- kappa_changes_with_ell=True
- epsilon_changes_with_ell=True

## Normalization rows
- lambda_sm_eff: accepted=False; expected_mass_dimension_4d=0; unknowns=['ell', 'A_action', 'Z_phi']
- kappa_gr_eff: accepted=False; expected_mass_dimension_4d=2; unknowns=['ell', 'A_action', 'Z_g']
- epsilon_mix_eff: accepted=False; expected_mass_dimension_4d=1; unknowns=['ell', 'A_action', 'Z_psi', 'Z_R']

## Missing normalization atoms
- canonical selector-free physical length/reference cell
- action-density unit convention for converting dimensionless sums into integral density
- cell volume / measure normalization for strict kernel moments M0..M3
- field and curvature normalization compatible with the same cell
- variational insertion check after fixing the above normalizations

## Decision
The audit finds a real scale-orbit/nonuniqueness witness: changing the positive reference length changes dimensionful kappa/epsilon while leaving the numeric moment map intact. P2689 and P2761 already block an unconditional unit/reference source, and no new canonical reference-cell/action-density theorem is exported.

## Recommendation
Do not promote the strict moment-map coefficients to physical Lagrangian couplings through an arbitrary unit choice.  The next honest proof-grade move should target exactly one remaining provenance atom after P2762: a sign-convention theorem for the moment-to-coupling map, but only if it is explicitly conditional on a future canonical reference-cell/action-density theorem; otherwise preserve the P2697-P2762 no-closure certificate.
