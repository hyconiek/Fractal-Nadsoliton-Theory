# P2415 S1365: phase/frequency affine-transport nonautomorphism certificate

Status: `PASS_PHASE_FREQUENCY_AFFINE_TRANSPORT_WITNESS_NO_SOURCE_NO_SELECTOR_CLOSURE`

## Result

P2415/S1365 packages the finite phase/frequency affine-coordinate transport and rejects Z12 reindexing or scalar replacement overreads.

## Finite facts

- Domain size: `12`.
- Max affine residual: `2.220446049250313e-16`.
- Z12 unit+offset checks: `48`.
- Best Z12 mismatch count: `2`.
- Phase bits: `[0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]`.
- Scalar replacement max residual: `0.9591793305520859`.

## Hard limits

- No strict dynamic source theorem for omega/phi is exported.
- No Z12 automorphism or scalar replacement of phase/frequency data is licensed.
- No orientation selector source or QW-2191 discharge follows.
- No full bridge completion, legacy role transfer, role-bearing L_total, or ToE closure follows.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'domain_has_twelve_nodes': True, 'affine_transport_residual_small': True, 'not_all_affine_coordinates_integer': True, 'checked_all_z12_unit_offsets': True, 'no_z12_reindex_match': True, 'scalar_phase_replacement_rejected': True, 'phase_bits_match_z2': True, 'gf2_unique_solution_inherited': True, 'no_phase_source_exported': True, 'no_selector_source_exported': True, 'phase_bridge_component_still_open': True, 'raw_identity_not_claimed': True, 'full_bridge_still_open': True, 'role_transfer_still_blocked': True, 'p2411_bridge_obligation_inherited': True, 'p2412_scope_separation_inherited': True, 'p2413_amplitude_witness_inherited': True, 'p2414_damping_nonabsorption_inherited': True, 'scratch_affine_inherited': True, 'scratch_nonautomorphism_inherited': True, 'fingerprint_stable': True}`
