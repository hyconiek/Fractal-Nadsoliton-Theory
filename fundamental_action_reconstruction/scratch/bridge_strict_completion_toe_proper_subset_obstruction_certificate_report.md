# ToE proper-subset obstruction certificate

Status: `PASS`

## Summary

- `all_127_proper_subsets_fail_toe`: `True`
- `nearest_miss_count`: `7`
- `six_atom_packages_fail_toe`: `True`
- `all_seven_atoms_required`: `True`
- `max_proper_true_targets`: `2`
- `toe_closure_claimed`: `False`

## Proper subset counts by weight

- weight `0`: `1` subsets
- weight `1`: `7` subsets
- weight `2`: `21` subsets
- weight `3`: `35` subsets
- weight `4`: `35` subsets
- weight `5`: `21` subsets
- weight `6`: `7` subsets

## Six-atom nearest misses

- missing `alpha_geo_electroweak_role_theorem` -> signature `1010`; role-transfer atom missing: electroweak role theorem is absent, so role transfer and ToE fail
- missing `beta_power_hierarchy_successor_theorem` -> signature `1010`; role-transfer atom missing: hierarchy successor theorem is absent, so role transfer and ToE fail
- missing `beta_tors_strict_role_theorem` -> signature `1010`; role-transfer atom missing: beta_tors strict successor theorem is absent, so role transfer and ToE fail
- missing `chi11_selector_source` -> signature `1000`; selector and role atom missing: QW-2191 selector, role transfer, and ToE fail
- missing `strict_damping_beta_eta_source` -> signature `0110`; bridge-source atom missing: damping beta/eta source is absent, so bridge and ToE fail
- missing `strict_dynamical_source_for_A_P_D` -> signature `0110`; bridge-source atom missing: APD source is absent, so bridge and ToE fail
- missing `strict_phase_frequency_source` -> signature `0110`; bridge-source atom missing: phase/frequency source is absent, so bridge and ToE fail

## Cross-checks

- `source_reports_present`: `True`
- `truth_table_count_pass`: `True`
- `proper_subset_count_pass`: `True`
- `full_set_only_toe_pass`: `True`
- `nearest_miss_count_pass`: `True`
- `nearest_miss_rows_fail_toe`: `True`
- `missing_atom_rows_cover_frontier`: `True`
- `minimal_cut_matches_truth_table`: `True`
- `target_lattice_full_signature_only`: `True`
- `readiness_certificate_inherited`: `True`
- `toe_audit_doc_mentions_joint_requirement`: `True`
- `hard_limits_preserved`: `True`

## Hard limits

- No proper subset is promoted to ToE closure.
- No new theorem atom is exported.
- No bridge theorem is claimed.
- No role-transfer theorem is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
