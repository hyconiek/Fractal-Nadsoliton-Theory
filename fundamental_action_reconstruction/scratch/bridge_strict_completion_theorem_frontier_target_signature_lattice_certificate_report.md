# Strict-completion theorem frontier target-signature lattice certificate

Status: `target-signature-lattice-enumerated-six-of-sixteen-signatures-reachable-no-closure`

## Result

The 128 open-atom assignments project to exactly 6 reachable target
signatures out of 16.  This records closure implications between the
targets without exporting any missing theorem atom.

## Summary

- `target_count`: `4`
- `all_target_signature_count`: `16`
- `reachable_target_signature_count`: `6`
- `unreachable_target_signature_count`: `10`
- `reachable_signatures`: `['0000', '0010', '0110', '1000', '1010', '1111']`
- `unreachable_signatures`: `['0001', '0011', '0100', '0101', '0111', '1001', '1011', '1100', '1101', '1110']`
- `counts_by_reachable_signature`: `{'0000': 56, '0010': 49, '0110': 7, '1000': 8, '1010': 7, '1111': 1}`
- `minimal_weights_by_reachable_signature`: `{'0000': 0, '0010': 1, '0110': 4, '1000': 3, '1010': 4, '1111': 7}`
- `role_implies_selector`: `True`
- `toe_implies_all_targets`: `True`
- `only_full_signature_has_toe_closure`: `True`
- `current_signature_all_false`: `True`
- `atom_influence_top_atom_inherited`: `True`
- `frontier_cut_open_atoms_match`: `True`
- `bridge_theorem_exported`: `False`
- `role_transfer_theorem_exported`: `False`
- `selector_closure_exported`: `False`
- `toe_closure_claimed`: `False`

## Signature rows

- `0000`: reachable=`True`, count=`56`, min_weight=`0`
- `0001`: reachable=`False`, count=`0`, min_weight=`None`
- `0010`: reachable=`True`, count=`49`, min_weight=`1`
- `0011`: reachable=`False`, count=`0`, min_weight=`None`
- `0100`: reachable=`False`, count=`0`, min_weight=`None`
- `0101`: reachable=`False`, count=`0`, min_weight=`None`
- `0110`: reachable=`True`, count=`7`, min_weight=`4`
- `0111`: reachable=`False`, count=`0`, min_weight=`None`
- `1000`: reachable=`True`, count=`8`, min_weight=`3`
- `1001`: reachable=`False`, count=`0`, min_weight=`None`
- `1010`: reachable=`True`, count=`7`, min_weight=`4`
- `1011`: reachable=`False`, count=`0`, min_weight=`None`
- `1100`: reachable=`False`, count=`0`, min_weight=`None`
- `1101`: reachable=`False`, count=`0`, min_weight=`None`
- `1110`: reachable=`False`, count=`0`, min_weight=`None`
- `1111`: reachable=`True`, count=`1`, min_weight=`7`

## Implications

- `role_transfer_theorem_level_closure => selector_qw2191_closure`: `True`
- `toe_closure => bridge_theorem_level_closure`: `True`
- `toe_closure => role_transfer_theorem_level_closure`: `True`
- `toe_closure => selector_qw2191_closure`: `True`

## Cross-checks

- `source_reports_present`: `True`
- `signature_partition_pass`: `True`
- `reachable_unreachable_counts_pass`: `True`
- `reachable_signature_counts_pass`: `True`
- `minimal_weights_pass`: `True`
- `implications_pass`: `True`
- `limits_preserved`: `True`

## Proof certificate

- `grep_step`: rg was used to avoid duplicating a target-signature or implication-lattice audit; none existed for the strict-completion theorem frontier.
- `projection_step`: Each of the 128 open-atom assignments is projected to a 4-bit target signature ordered as bridge, role-transfer, selector, ToE.
- `reachability_step`: Exactly 6 of 16 target signatures are reachable: 0000, 0010, 0110, 1000, 1010, and 1111.
- `count_step`: The reachable signature counts are 56, 49, 7, 8, 7, and 1 respectively, summing to 128.
- `implication_step`: The lattice proves role-transfer closure implies selector closure, and ToE closure implies bridge, role-transfer, and selector closure.
- `minimal_step`: Minimal true-atom weights by reachable signature are 0, 1, 4, 3, 4, and 7 respectively; ToE is reachable only at the all-open-atom frontier cut.
- `scope_step`: This is a finite target-lattice audit only; it exports no theorem atom and proves no bridge, role-transfer, selector, or ToE closure.

## Hard limits

- No reachable target signature is promoted to the current theory state.
- No missing theorem atom is exported.
- No bridge theorem, role-transfer theorem, selector closure, or ToE closure is claimed.
- No QW-2191 selector discharge is claimed.
