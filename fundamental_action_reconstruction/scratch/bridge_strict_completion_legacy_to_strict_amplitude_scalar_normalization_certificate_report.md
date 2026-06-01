# Legacy-to-strict amplitude scalar-normalization certificate

Status: `legacy-alpha-geo-factors-out-as-positive-global-scalar__bridge-row-witness-not-role-transfer`

## Result

The visible legacy scalar `alpha_geo=4 ln(2)` factors out of `K_legacy_ont`
as a positive global amplitude scalar on the audited finite domain.  This is
a scalar-normalization witness only, not a physical-role transfer theorem.

## Summary

- `legacy_alpha_geo_visible`: `True`
- `strict_alpha_geo_source_loaded`: `True`
- `finite_domain_nodes`: `[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]`
- `legacy_shape_nonzero_on_domain`: `True`
- `legacy_denominator_positive_on_domain`: `True`
- `cos_zero_congruence_has_no_domain_solution`: `True`
- `alpha_geo_positive`: `True`
- `alpha_inverse_normalization_residual_zero_formally`: `True`
- `alpha_inverse_normalization_residual_max_abs_float`: `1.1102230246251565e-16`
- `ratio_K_over_shape_constant_alpha_geo_formally`: `True`
- `ratio_K_over_shape_max_float_deviation`: `4.440892098500626e-16`
- `sign_pattern_preserved_by_positive_alpha`: `True`
- `necessity_exact_apd_inherited`: `True`
- `component_gap_amplitude_row_present`: `True`
- `scalar_normalization_witness_exported`: `True`
- `full_strict_A_d_derivation_exported`: `False`
- `legacy_role_transfer_allowed`: `False`
- `raw_kernel_identity_claimed`: `False`
- `full_bridge_theorem_exported`: `False`

## Cross-checks

- `guardrail_requires_completion_map`: `True`
- `component_gap_amplitude_row_open_but_finite`: `True`
- `alpha_source_strict_derived_value_loaded`: `True`
- `scalar_normalization_formally_exact`: `True`
- `finite_domain_safe_for_ratio`: `True`
- `signs_preserved`: `True`
- `no_role_transfer_or_bridge_claim`: `True`

## Proof certificate

- `grep_step`: rg was used to distinguish this scalar-normalization bridge witness from older amplitude nonabsorption/role-transfer packets.
- `factorization_step`: For every audited d, K_legacy_ont(d)=alpha_geo*L_shape(d) by direct algebra; applying alpha_geo^{-1} gives exactly L_shape(d).
- `nonzero_step`: On d=0..11, cos(pi*d/4+pi/6)=cos((3d+2)pi/12) cannot vanish because 3d==4 mod 12 has no solution; the denominator 1+0.01*d is positive.
- `positivity_step`: alpha_geo=4 ln(2)>0, so scalar normalization preserves the legacy sign pattern and cannot supply the strict GF(2) phase flips.
- `strict_source_step`: The repo already exports alpha_geo_strict_derived_v1=4 ln(2), but this certificate uses it only as scalar provenance and not as a selector or physical-role transfer theorem.
- `bridge_meaning_step`: The amplitude row now has an explicit scalar-normalization witness; the remaining bridge still needs strict A(d) source, phase/frequency transport source, nonlinear damping/compression source, selector source, and role-transfer audit.
- `theoretical_limit`: This is not a full A_abs absorption theorem, not a full legacy->strict bridge, and not a transfer of sin^2(theta_W), alpha_EM, beta^N, or beta_tors->chi_11 roles.

## Hard limits

- No raw identity K_legacy_ont == K_strict_gate is claimed.
- No full strict A(d) dynamical derivation is exported.
- No legacy alpha_geo physical-role transfer is licensed.
- No beta_tors -> chi_11 theorem is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
