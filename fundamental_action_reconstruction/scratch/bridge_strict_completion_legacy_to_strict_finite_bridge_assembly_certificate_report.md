# Legacy-to-strict finite bridge assembly certificate

Status: `finite-apd-bridge-assembly-reconstructs-strict-kernel-on-z12__source-selector-role-transfer-open`

## Result

The separate amplitude, phase/frequency, and damping/compression witnesses
assemble into one finite kernel-comparison map `Q_assembly=A*P*D` on Z12.
The assembled map reconstructs `K_strict_gate` from `K_legacy_ont` on the
audited finite domain, but remains below source/selector/role-transfer closure.

## Summary

- `domain_size`: `12`
- `finite_kernel_comparison_witness_exported`: `True`
- `t99_positive_bridge_closure_target_detected`: `True`
- `comparison_scope_only`: `True`
- `amplitude_scalar_normalization_inherited`: `True`
- `phase_affine_transport_inherited`: `True`
- `damping_beta_eta_identifiability_inherited`: `True`
- `finite_diagonal_completion_map_inherited`: `True`
- `max_abs_reconstruction_residual`: `2.220446049250313e-16`
- `max_abs_assembled_q_minus_diagonal_q`: `1.1102230246251565e-16`
- `max_abs_assembled_q_minus_necessity_factor_product`: `0.0`
- `max_abs_log_additive_residual`: `8.881784197001252e-16`
- `assembled_map_reconstructs_strict_kernel`: `True`
- `assembled_map_matches_finite_diagonal_certificate`: `True`
- `assembled_map_matches_necessity_apd_product`: `True`
- `log_abs_decomposition_additive`: `True`
- `phase_sign_bits_match_diagonal_sign_bits`: `True`
- `phase_sign_bits_match_prior_affine_report`: `True`
- `component_gap_still_not_full_bridge`: `True`
- `guardrail_role_transfer_required_after_full_bridge`: `True`
- `strict_dynamic_source_exported`: `False`
- `legacy_role_transfer_exported`: `False`
- `selector_source_exported`: `False`
- `full_bridge_theorem_exported`: `False`
- `toe_closure_claimed`: `False`

## Cross-checks

- `source_reports_present`: `True`
- `t99_target_scope_detected`: `True`
- `component_certificates_inherited`: `True`
- `assembled_identity_exact_on_finite_domain`: `True`
- `sign_and_log_decomposition_coherent`: `True`
- `guardrail_limits_preserved`: `True`

## Proof certificate

- `grep_step`: rg was used before adding this probe to avoid duplicating the separate amplitude, phase-affine, damping-parameter, and diagonal certificates.
- `assembly_step`: For every audited d, Q_assembly(d)=alpha_geo^{-1} * cos(theta_S(d))/cos(theta_L(d)) * (1+beta_tors*d)/(1+beta*d^eta).
- `identity_step`: Substituting K_legacy_ont into K_legacy_ont(d)*Q_assembly(d) cancels alpha_geo, cos(theta_L), and the legacy denominator, leaving K_strict_gate(d).
- `diagonal_step`: The assembled Q matches the previously certified unique diagonal Q=diag(K_strict/K_legacy) on all 12 audited nodes.
- `log_sign_step`: Since A and D are positive, sign(Q_assembly) is carried by P; log|Q|=log(A)+log|P|+log(D) numerically closes on every audited node.
- `scope_step`: This is a finite T99-style kernel-comparison witness only; it does not export a strict dynamical source, selector source, role-transfer theorem, QW-2191 discharge, or ToE closure.

## Hard limits

- No unqualified raw identity K_legacy_ont == K_strict_gate is claimed.
- No strict dynamical derivation of A/P/D, omega/phi, beta/eta, or chi_11 is exported.
- No beta_tors -> beta/eta or beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer is licensed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
