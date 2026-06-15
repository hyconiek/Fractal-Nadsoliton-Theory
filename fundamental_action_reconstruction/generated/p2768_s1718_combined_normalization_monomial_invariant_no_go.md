# P2768/S1718 combined-normalization monomial-invariant no-go

Status: `P2768_COMBINED_NORMALIZATION_MONOMIAL_INVARIANT_NO_GO_NO_CLOSURE`

## Weight matrix
- coefficient_order=['lambda_sm_eff', 'kappa_gr_eff', 'epsilon_mix_eff']
- gauge_order=['ell_reference_length', 'Z_phi_scalar_field', 'Z_R_curvature']
- weight_matrix=[[0, -2, -1], [-4, 0, -2], [0, -1, -1]]
- determinant=4
- full_rank_over_Q=True

## Brute-force invariant scan
- tested_nonzero_exponent_vectors=728
- nontrivial_invariant_count=0

## Ratio candidates
- epsilon_squared_over_lambda_kappa: weight=(0, 0, -1); invariant=False; distinct_values=3
- epsilon_squared_over_lambda_kappa_squared: weight=(2, 0, 0); invariant=False; distinct_values=4
- kappa_over_epsilon_squared: weight=(0, 4, 1); invariant=False; distinct_values=4
- p1562_R1_R2_R3_moment_ratios_not_coupling_normalization: weight=(0, 0, 0); invariant=True; distinct_values=1

## Decision
The combined length/field/curvature normalization action has full-rank weight matrix on the three P1562 coefficients.  Therefore no nontrivial monomial ratio of lambda_sm_eff, kappa_gr_eff, and epsilon_mix_eff is invariant under the combined open normalizations; this blocks a ratio-invariant rescue but does not export a canonical normalization theorem.

## Recommendation
Do not try to rescue physical-coupling provenance by taking a monomial ratio of lambda_sm_eff, kappa_gr_eff, and epsilon_mix_eff.  P2768 proves the combined open normalization action has no nontrivial monomial invariant.  The next honest step must either supply an actual canonical normalization theorem/source, or introduce a genuinely non-monomial invariant with an explicit invariance proof and then run the same bounded acceptance test; otherwise preserve the P2697-P2768 no-closure certificate.
