# P2405 S1355: nadsoliton information-ontology projection certificate

Status: `PASS_INFORMATION_ONTOLOGY_RETYPE_NO_UNDERLAYER_NO_ROLE_EXPORT`

## Result

P2405/S1355 retypes P2404 strict additions as internal information constraints of the one nadsoliton, not as a lower information layer or direct physical-role exports.

## Exact ontology guard

- Guard true masks: `[31]`.
- Proper-subset fail count: `31`.
- Guard ANF: `nadsoliton_is_sole_primordial_information * no_separate_information_layer_under_nadsoliton * strict_additions_are_internal_information_constraints * physics_roles_are_downstream_projections * observer_is_downstream_readout_not_source`.
- Guard ANF degree: `5`.

## Poset / projection facts

- Unique information root: `True`.
- All nodes reachable from information root: `True`.
- Topological order count: `1`.
- P2404 strict-additions-only physical role lanes ready: `[]`.
- P2404 L_total dependency degree: `7`.

## Hard limits

- No separate information substrate below the nadsoliton is introduced.
- No strict addition is retyped as an immediate physical role theorem.
- No observer/readout source closes QW-2191 or supplies a strict selector theorem.
- No L_total or ToE promotion follows from ontology typing alone.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'ax9_single_primitive_inherited': True, 'ax9_blocks_separate_information_layer': True, 'ontology_lattice_has_32_rows': True, 'only_full_ontology_guard_mask_passes': True, 'ontology_guard_degree_five': True, 'unique_nadsoliton_information_root': True, 'all_nodes_reachable_from_information_root': True, 'p2404_strict_additions_do_not_license_physical_roles': True, 'p2404_degree_seven_dependency_inherited': True, 'n50_internal_route_warning_inherited': True, 's2_priority_order_inherited': True, 'fingerprint_stable': True}`
