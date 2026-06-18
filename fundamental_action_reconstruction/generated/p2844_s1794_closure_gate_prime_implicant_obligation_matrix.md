# P2844/S1794 closure-gate prime-implicant obligation matrix

Status: `P2844_CLOSURE_GATE_PRIME_IMPLICANT_OBLIGATION_MATRIX_NO_CLOSURE`

## Current evidence atoms
`finite_witness_available`, `kernel_split_preserved`, `typed_source_domain`

## Target obstruction matrix
### strict_kernel_completion_bridge
- currently_closed=False
- missing_count=4
- missing_atoms=['bridge_amplitude_map', 'bridge_damping_compression_map', 'bridge_phase_frequency_map', 'bridge_selector_source']
- minimal_unlock_sets=[['bridge_amplitude_map', 'bridge_damping_compression_map', 'bridge_phase_frequency_map', 'bridge_selector_source']]
### legacy_role_transfer
- currently_closed=False
- missing_count=5
- missing_atoms=['bridge_amplitude_map', 'bridge_damping_compression_map', 'bridge_phase_frequency_map', 'bridge_selector_source', 'role_transfer_theorem']
- minimal_unlock_sets=[['bridge_amplitude_map', 'bridge_damping_compression_map', 'bridge_phase_frequency_map', 'bridge_selector_source', 'role_transfer_theorem']]
### typed_ltotal_source_coupling
- currently_closed=False
- missing_count=7
- missing_atoms=['coupling_coefficient_rule', 'locality_covariance', 'localization_pullback', 'nonproxy_ltotal_term', 'target_independent_units', 'typed_target_codomain', 'variational_chain_rule']
- minimal_unlock_sets=[['coupling_coefficient_rule', 'locality_covariance', 'localization_pullback', 'nonproxy_ltotal_term', 'target_independent_units', 'typed_target_codomain', 'variational_chain_rule']]
### eom_closure
- currently_closed=False
- missing_count=8
- missing_atoms=['coupling_coefficient_rule', 'eom_boundary_rules', 'localization_pullback', 'nonproxy_eom_residual_zero', 'nonproxy_ltotal_term', 'target_independent_units', 'typed_target_codomain', 'variational_chain_rule']
- minimal_unlock_sets=[['coupling_coefficient_rule', 'eom_boundary_rules', 'localization_pullback', 'nonproxy_eom_residual_zero', 'nonproxy_ltotal_term', 'target_independent_units', 'typed_target_codomain', 'variational_chain_rule']]
### hamiltonian_closure
- currently_closed=False
- missing_count=7
- missing_atoms=['eom_boundary_rules', 'hamiltonian_boundedness', 'hamiltonian_recovers_eom', 'legendre_regular_or_constraint_split', 'nonproxy_eom_residual_zero', 'nonproxy_ltotal_term', 'variational_chain_rule']
- minimal_unlock_sets=[['eom_boundary_rules', 'hamiltonian_boundedness', 'hamiltonian_recovers_eom', 'legendre_regular_or_constraint_split', 'nonproxy_eom_residual_zero', 'nonproxy_ltotal_term', 'variational_chain_rule']]
### toe_style_promotion
- currently_closed=False
- missing_count=18
- missing_atoms=['bridge_amplitude_map', 'bridge_damping_compression_map', 'bridge_phase_frequency_map', 'bridge_selector_source', 'coupling_coefficient_rule', 'eom_boundary_rules', 'hamiltonian_boundedness', 'hamiltonian_recovers_eom', 'legendre_regular_or_constraint_split', 'locality_covariance', 'localization_pullback', 'nonproxy_eom_residual_zero', 'nonproxy_ltotal_term', 'role_transfer_theorem', 'selector_qw2191_discharged', 'target_independent_units', 'typed_target_codomain', 'variational_chain_rule']
- minimal_unlock_sets=[['bridge_amplitude_map', 'bridge_damping_compression_map', 'bridge_phase_frequency_map', 'bridge_selector_source', 'coupling_coefficient_rule', 'eom_boundary_rules', 'hamiltonian_boundedness', 'hamiltonian_recovers_eom', 'legendre_regular_or_constraint_split', 'locality_covariance', 'localization_pullback', 'nonproxy_eom_residual_zero', 'nonproxy_ltotal_term', 'role_transfer_theorem', 'selector_qw2191_discharged', 'target_independent_units', 'typed_target_codomain', 'variational_chain_rule']]

## Candidate next-move scores
- unit_bearing_typed_source_coupling_map: admissible=True, directly_closes=[], missing_after={'strict_kernel_completion_bridge': 4, 'legacy_role_transfer': 5, 'typed_ltotal_source_coupling': 2, 'eom_closure': 4, 'hamiltonian_closure': 7, 'toe_style_promotion': 13}
- variational_chain_rule_theorem: admissible=False, directly_closes=[], missing_after={'strict_kernel_completion_bridge': 4, 'legacy_role_transfer': 5, 'typed_ltotal_source_coupling': 6, 'eom_closure': 6, 'hamiltonian_closure': 5, 'toe_style_promotion': 16}
- selector_source_provider: admissible=True, directly_closes=[], missing_after={'strict_kernel_completion_bridge': 3, 'legacy_role_transfer': 4, 'typed_ltotal_source_coupling': 7, 'eom_closure': 8, 'hamiltonian_closure': 7, 'toe_style_promotion': 16}
- single_kernel_bridge_atom_amplitude: admissible=True, directly_closes=[], missing_after={'strict_kernel_completion_bridge': 3, 'legacy_role_transfer': 4, 'typed_ltotal_source_coupling': 7, 'eom_closure': 8, 'hamiltonian_closure': 7, 'toe_style_promotion': 17}
- single_kernel_bridge_atom_damping: admissible=True, directly_closes=[], missing_after={'strict_kernel_completion_bridge': 3, 'legacy_role_transfer': 4, 'typed_ltotal_source_coupling': 7, 'eom_closure': 8, 'hamiltonian_closure': 7, 'toe_style_promotion': 17}
- hamiltonian_legendre_only: admissible=False, directly_closes=[], missing_after={'strict_kernel_completion_bridge': 4, 'legacy_role_transfer': 5, 'typed_ltotal_source_coupling': 7, 'eom_closure': 8, 'hamiltonian_closure': 6, 'toe_style_promotion': 17}

## Boundary
P2844 exhaustively computes the Boolean closure-obligation matrix.  The only currently true atoms are finite witness availability, kernel-split preservation, and a finite source-domain placeholder.  Every promoted target has nonempty missing prime-implicant atoms, so no bridge, role transfer, L_total, EOM, Hamiltonian, selector, or ToE closure is exported.

## Recommendation
Attack one minimal high-leverage atom bundle rather than replaying SWOT or graph separation: first choice is a unit-bearing typed source/coupling map exporting target codomain, units, localization/pullback, locality/covariance, and coupling coefficient; second choice is exactly one kernel bridge atom, preferably damping/compression or amplitude, with an explicit source premise.  A Hamiltonian-only Legendre move is not admissible before nonproxy L_total/EOM closure.
