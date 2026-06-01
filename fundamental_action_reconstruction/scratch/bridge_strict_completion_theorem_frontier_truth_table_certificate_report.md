# Strict-completion theorem frontier truth-table certificate

Status: `all-2pow7-frontier-assignments-enumerated-current-assignment-closes-no-target`

## Result

All 128 assignments of the seven open theorem atoms are enumerated.
The current all-false assignment closes no target, and ToE closure requires
all seven atoms.

## Summary

- `open_atom_count`: `7`
- `truth_assignment_count`: `128`
- `current_assignment_all_false`: `True`
- `current_targets_all_false`: `True`
- `bridge_satisfying_assignment_count`: `16`
- `role_satisfying_assignment_count`: `8`
- `selector_satisfying_assignment_count`: `64`
- `toe_satisfying_assignment_count`: `1`
- `bridge_minimal_set_size`: `3`
- `role_minimal_set_size`: `4`
- `selector_minimal_set_size`: `1`
- `toe_minimal_set_size`: `7`
- `toe_minimal_set_equals_frontier_cut`: `True`
- `role_lattice_minimal_set_inherited`: `True`
- `component_gap_sources_still_missing`: `True`
- `anchor_selector_source_still_open`: `True`
- `bridge_theorem_exported`: `False`
- `role_transfer_theorem_exported`: `False`
- `selector_closure_exported`: `False`
- `toe_closure_claimed`: `False`

## Cross-checks

- `source_reports_present`: `True`
- `truth_table_size_pass`: `True`
- `current_assignment_no_closure`: `True`
- `minimal_sizes_match_definitions`: `True`
- `toe_minimal_set_matches_frontier_cut`: `True`
- `limits_preserved`: `True`

## Proof certificate

- `grep_step`: rg was used to avoid duplicating the theorem-frontier cut; this report adds the exhaustive Boolean truth table over its seven open atoms.
- `enumeration_step`: All 2^7=128 truth assignments of the open theorem atoms are enumerated exactly.
- `target_step`: Bridge closure requires three strict-source atoms; role-transfer closure requires four role/selector atoms; selector closure requires chi11; ToE requires all target predicates.
- `current_step`: The current assignment has all seven atoms false, so bridge, role-transfer, selector, and ToE targets are all false.
- `minimal_step`: The minimal true-atom sizes are bridge=3, role=4, selector=1, and ToE=7, matching the frontier-cut report.
- `scope_step`: This is a Boolean readiness certificate only; it exports no theorem atom, no QW-2191 discharge, and no ToE closure.

## Hard limits

- No truth-table assignment is promoted to the current theory state.
- No missing theorem atom is exported.
- No bridge theorem, role-transfer theorem, selector closure, or ToE closure is claimed.
- No QW-2191 selector discharge is claimed.
