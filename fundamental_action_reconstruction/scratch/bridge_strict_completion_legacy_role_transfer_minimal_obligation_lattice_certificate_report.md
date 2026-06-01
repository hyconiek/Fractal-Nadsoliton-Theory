# Legacy role-transfer minimal obligation lattice certificate

Status: `minimal-role-transfer-obligation-lattice-enumerated-all-theorem-atoms-missing`

## Result

All subsets of the missing role-transfer theorem atoms are enumerated.
The unique minimal set covering all audited legacy roles contains all four
atoms, and none is exported by this certificate.

## Summary

- `role_count`: `4`
- `theorem_atom_count`: `4`
- `total_subset_count_checked`: `16`
- `all_pre_audit_roles_loaded`: `True`
- `all_roles_blocked_in_pre_audit`: `True`
- `all_atoms_missing_now`: `True`
- `global_minimal_obligation_count`: `1`
- `global_minimal_obligation_size`: `4`
- `global_minimal_obligation_sets`: `[['alpha_geo_electroweak_role_theorem', 'beta_tors_strict_role_theorem', 'beta_power_hierarchy_successor_theorem', 'chi11_selector_source_theorem']]`
- `beta_tors_atom_is_shared_by_three_roles`: `True`
- `symbolic_bridge_still_no_role_transfer`: `True`
- `component_gap_still_blocks_role_transfer`: `True`
- `guardrail_requires_audit`: `True`
- `role_transfer_theorem_exported`: `False`
- `q2191_discharged`: `False`
- `toe_closure_claimed`: `False`

## Role lattice rows

- `legacy_weak_mixing_angle` requires `['alpha_geo_electroweak_role_theorem']`
- `legacy_inverse_alpha_em` requires `['alpha_geo_electroweak_role_theorem', 'beta_tors_strict_role_theorem']`
- `legacy_beta_power_gravity_hierarchy` requires `['beta_power_hierarchy_successor_theorem', 'beta_tors_strict_role_theorem']`
- `legacy_torsion_to_chi11_orientation` requires `['beta_tors_strict_role_theorem', 'chi11_selector_source_theorem']`

## Cross-checks

- `source_reports_present`: `True`
- `pre_audit_roles_and_blocks_inherited`: `True`
- `global_minimal_set_is_all_atoms`: `True`
- `beta_tors_shared_obligation_detected`: `True`
- `all_atoms_missing_and_limits_preserved`: `True`

## Proof certificate

- `grep_step`: rg was used to distinguish this minimal-obligation lattice from the earlier pre-audit and obstruction notes.
- `enumeration_step`: All 16 subsets of the four theorem atoms were enumerated.
- `role_step`: Each role is satisfied only by subsets containing its required theorem atoms; with zero atoms exported now, every transfer remains blocked.
- `global_step`: To cover all four audited roles, the unique minimal obligation set is the full four-atom set: alpha_geo role, beta_tors role, beta-power successor, and chi11 selector/source.
- `shared_beta_step`: The beta_tors strict-role theorem is shared by alpha_EM, beta^N hierarchy, and beta_tors->chi11 obligations, so it is the widest missing role-transfer bottleneck.
- `scope_step`: This lattice is a planning certificate only; it exports no role-transfer theorem, no QW-2191 discharge, and no ToE closure.

## Hard limits

- No theorem atom in the role-transfer obligation lattice is exported here.
- No legacy physical-role claim is transferred onto K_strict_gate.
- No alpha_geo, beta_tors, beta^N, or chi11 role theorem is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
