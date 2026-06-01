# Legacy-to-strict phase/frequency affine transport certificate

Status: `continuous-affine-phase-transport-exact__not-z12-automorphism__phase-source-open`

## Result

The continuous affine coordinate `x(d)` exactly transports the legacy phase
argument to the strict phase argument, but it is not a `Z12` automorphism
or scalar replacement.  The selector/source question remains open.

## Summary

- `domain_size`: `12`
- `affine_slope_omega_s_over_omega_l`: `0.23650424543455648`
- `affine_intercept`: `-0.4597652406472027`
- `continuous_affine_phase_transport_exact`: `True`
- `max_abs_affine_transport_residual`: `2.220446049250313e-16`
- `affine_map_is_not_z12_automorphism`: `True`
- `affine_coordinates_not_all_integers`: `True`
- `minimum_distance_to_integer_affine_coordinate`: `0.01324325022191032`
- `maximum_distance_to_integer_affine_coordinate`: `0.4862517410910232`
- `z12_unit_offset_automorphism_count_checked`: `48`
- `no_z12_unit_offset_reindex_matches_strict_sign_pattern`: `True`
- `best_z12_unit_offset_mismatch_count`: `2`
- `phase_factor_bits`: `[0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]`
- `phase_factor_bits_match_z2_node_bits`: `True`
- `gf2_solution_inherited_unique`: `True`
- `scalar_phase_replacement_fails`: `True`
- `scalar_phase_best_fit_max_abs_residual`: `0.9591793305520859`
- `scalar_phase_best_fit_l2_residual`: `2.1121772977850153`
- `finite_diagonal_completion_map_inherited`: `True`
- `component_gap_phase_row_still_source_open`: `True`
- `strict_phase_frequency_source_exported`: `False`
- `orientation_selector_source_exported`: `False`
- `raw_kernel_identity_claimed`: `False`
- `full_bridge_theorem_exported`: `False`

## Cross-checks

- `guardrail_restores_legacy_as_intermediate`: `True`
- `continuous_affine_transport_exact`: `True`
- `not_discrete_z12_automorphism`: `True`
- `phase_bits_match_gf2_chain`: `True`
- `scalar_phase_replacement_rejected`: `True`
- `no_source_selector_or_bridge_claim`: `True`

## Proof certificate

- `grep_step`: rg was used to separate this affine phase-frequency transport certificate from existing pointwise phase-factor and GF(2) sign reports.
- `affine_step`: By construction theta_L(x(d))=omega_L*((omega_S*d+phi_S-phi_L)/omega_L)+phi_L=theta_S(d), so the continuous phase-coordinate transport is exact.
- `non_automorphism_step`: The affine slope omega_S/omega_L is not a Z12 unit and the sampled x(d) are not all integers; exhaustive unit+offset checks over Aut(Z12) find no reindexing that reproduces the strict cosine sign pattern.
- `phase_factor_step`: The pointwise factor P(d)=cos(theta_S(d))/cos(theta_L(d)) has GF(2) sign bits matching the existing Z2/GF(2) phase-sign chain.
- `non_scalar_step`: The best scalar fit from the legacy cosine carrier to the strict cosine carrier has nonzero residual, so phase/frequency transport is not scalar normalization.
- `theoretical_limit`: This is a finite/affine phase-frequency bridge witness, not a strict derivation of omega/phi, not an orientation selector source, not beta_tors->chi_11, and not ToE closure.

## Hard limits

- No raw identity K_legacy_ont == K_strict_gate is claimed.
- No strict dynamical source for omega/phi or phase transport is exported.
- No orientation/selector source is exported.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer is licensed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
