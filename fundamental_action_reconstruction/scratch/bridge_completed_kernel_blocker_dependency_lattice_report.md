# Completed kernel blocker dependency lattice

Status: `minimal-premise-lattice-computed-for-completed-strict-kernel-bridge-frontier`

## Finite lattice

- Premises: `13`
- Enumerated subsets: `8192`
- Exported certificates: `6`
- Open blockers: `7`

## Current frontier summary

- `current_live_kernel`: K_strict_gate is the live full kernel; legacy is the historical/nadsoliton-characteristic carrier.
- `already_certified`: ['diagrams_legacy_carrier', 'exact_completion_factorization', 'local_puiseux_match', 'eta_plus_one_puiseux_match', 'measure_transport_identity', 'monotone_flow_output_matching']
- `remaining_for_theorem_level_bridge`: ['strict_transport_derivation', 'global_z12_map_derivation', 'orientation_chi11_source', 'chi11_uniqueness', 'reynolds_obstruction_escape', 'role_transfer_theorem']
- `main_one_bit_blocker`: orientation_chi11_source plus chi11_uniqueness and reynolds_obstruction_escape
- `transport_blocker`: strict_transport_derivation plus global_z12_map_derivation
- `role_and_selector_blocker`: role_transfer_theorem remains below any QW-2191 selector discharge

## Outcome antichains

### finite_exact_completion_certificate
- Minimal premise count: `1`
- Currently realized: `True`
- Missing blockers for first set: `[]`
- Minimal sets: `[['exact_completion_factorization']]`

### guarded_completed_kernel_candidate
- Minimal premise count: `2`
- Currently realized: `True`
- Missing blockers for first set: `[]`
- Minimal sets: `[['diagrams_legacy_carrier', 'exact_completion_factorization']]`

### local_two_channel_completion_candidate
- Minimal premise count: `1`
- Currently realized: `True`
- Missing blockers for first set: `[]`
- Minimal sets: `[['local_puiseux_match']]`

### eta_plus_one_local_completion_candidate
- Minimal premise count: `2`
- Currently realized: `True`
- Missing blockers for first set: `[]`
- Minimal sets: `[['local_puiseux_match', 'eta_plus_one_puiseux_match']]`

### measure_transport_bookkeeping_candidate
- Minimal premise count: `1`
- Currently realized: `True`
- Missing blockers for first set: `[]`
- Minimal sets: `[['measure_transport_identity']]`

### monotone_output_matching_flow_candidate
- Minimal premise count: `1`
- Currently realized: `True`
- Missing blockers for first set: `[]`
- Minimal sets: `[['monotone_flow_output_matching']]`

### theorem_level_completed_legacy_to_strict_bridge
- Minimal premise count: `8`
- Currently realized: `False`
- Missing blockers for first set: `['strict_transport_derivation', 'global_z12_map_derivation', 'orientation_chi11_source', 'chi11_uniqueness', 'reynolds_obstruction_escape', 'role_transfer_theorem']`
- Minimal sets: `[['diagrams_legacy_carrier', 'exact_completion_factorization', 'strict_transport_derivation', 'global_z12_map_derivation', 'orientation_chi11_source', 'chi11_uniqueness', 'reynolds_obstruction_escape', 'role_transfer_theorem']]`

### selector_closed_completed_kernel_toe_step
- Minimal premise count: `9`
- Currently realized: `False`
- Missing blockers for first set: `['strict_transport_derivation', 'global_z12_map_derivation', 'orientation_chi11_source', 'chi11_uniqueness', 'reynolds_obstruction_escape', 'role_transfer_theorem', 'qw2191_selector_discharge']`
- Minimal sets: `[['diagrams_legacy_carrier', 'exact_completion_factorization', 'strict_transport_derivation', 'global_z12_map_derivation', 'orientation_chi11_source', 'chi11_uniqueness', 'reynolds_obstruction_escape', 'role_transfer_theorem', 'qw2191_selector_discharge']]`

### strict_full_aut_internal_chi11_source
- Minimal premise count: `None`
- Currently realized: `False`
- Missing blockers for first set: `None`
- Minimal sets: `[]`

## Proof certificate

- `finite_domain`: The dependency lattice enumerates all 2^13=8192 subsets of named bridge premises.
- `separation_result`: Exact completion factorization, local Puiseux matching, measure transport bookkeeping, and monotone output-flow candidates are already certificates, but they are not strict derivations.
- `bridge_antichain`: The theorem-level completed legacy-to-strict bridge requires strict transport derivation, global Z12 map derivation, chi_11 source/uniqueness/Reynolds escape, and role-transfer theorem in addition to the existing carrier/factorization certificates.
- `full_aut_limit`: The strict_full_aut_internal_chi11_source outcome has no minimal premise set in this lattice because current Reynolds evidence keeps full-Aut invariant data in the annihilating/orthogonal sector.

## Hard limits

- No claim that the legacy kernel is the current live kernel; K_strict_gate is the current full form.
- No unqualified identity K_legacy_ont == K_strict_gate is asserted.
- No beta_tors -> chi_11 theorem is asserted.
- No strict derivation of the completion factors, transport ODE, or global Z12 map is claimed.
- No legacy physical-role transfer onto K_strict_gate is used without an explicit bridge theorem.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
