# Legacy physical-role transfer pre-audit certificate

Status: `legacy-physical-role-claims-audited-all-blocked-until-separate-role-transfer-theorem`

## Result

The listed legacy physical-role claims are dependency-audited and all remain
blocked until a separate role-transfer theorem is supplied.  The symbolic
bridge/cancellation certificates are kernel-comparison evidence only.

## Summary

- `role_claim_count`: `4`
- `dependency_matrix_columns`: `['depends_alpha_geo', 'depends_beta_tors', 'depends_beta_power_hierarchy', 'depends_chi11_or_orientation', 'strict_successor_theorem_available_now', 'unchanged_transfer_allowed_now', 'modified_transfer_allowed_now', 'rejected_by_current_bridge_data']`
- `dependency_matrix_rank_gf2`: `4`
- `all_roles_blocked_now`: `True`
- `roles_transferred_now`: `0`
- `s2_lists_required_role_transfer_claims`: `True`
- `t15_keeps_role_transfer_separate`: `True`
- `guardrail_requires_role_transfer_audit`: `True`
- `symbolic_bridge_is_formula_only`: `True`
- `component_gap_blocks_role_transfer`: `True`
- `alpha_role_source_missing`: `True`
- `beta_tors_role_source_missing`: `True`
- `chi11_source_missing`: `True`
- `strict_dynamic_source_exported`: `False`
- `role_transfer_theorem_exported`: `False`
- `q2191_discharged`: `False`
- `toe_closure_claimed`: `False`

## Role rows

- `legacy_weak_mixing_angle`: `blocked_not_transferred` — The bridge cancels alpha_geo as a kernel amplitude scalar; no theorem says the legacy electroweak role survives that normalization.
- `legacy_inverse_alpha_em`: `blocked_not_transferred` — The current bridge separates legacy linear beta_tors from strict nonlinear beta*d^eta and exports no beta_tors role theorem.
- `legacy_beta_power_gravity_hierarchy`: `blocked_not_transferred` — No report maps the legacy beta-power hierarchy through the strict beta=1, eta=9/5 nonlinear compression data.
- `legacy_torsion_to_chi11_orientation`: `blocked_not_transferred` — GF(2)/H1 reports locate the bit but do not source it, and the anchor/H1 audit keeps selector source open.

## Cross-checks

- `source_reports_present`: `True`
- `s2_and_t15_require_separate_role_transfer`: `True`
- `all_roles_blocked_with_nonzero_dependency_rank`: `True`
- `component_and_symbolic_reports_do_not_transfer_roles`: `True`
- `alpha_beta_chi11_sources_missing`: `True`
- `closure_limits_preserved`: `True`

## Proof certificate

- `grep_step`: rg was used to locate existing role-transfer guardrails and avoid treating the symbolic bridge as a physical-role theorem.
- `matrix_step`: The four audited legacy roles give a GF(2) dependency matrix of rank 4 over columns ['depends_alpha_geo', 'depends_beta_tors', 'depends_beta_power_hierarchy', 'depends_chi11_or_orientation', 'strict_successor_theorem_available_now', 'unchanged_transfer_allowed_now', 'modified_transfer_allowed_now', 'rejected_by_current_bridge_data'].
- `alpha_step`: sin^2(theta_W)=alpha_geo/12 is blocked because alpha_geo is only normalized as a kernel scalar, with no strict physical-role theorem.
- `beta_step`: alpha_EM^-1 and beta^N hierarchy claims are blocked because beta_tors is not mapped to a strict physical role through beta=1, eta=9/5 compression.
- `chi11_step`: beta_tors -> chi_11 is blocked because the anchor/H1 certificates locate the bit but keep its selector/source theorem open.
- `scope_step`: This is a pre-audit: it classifies transfer obligations but exports zero role-transfer permissions and no QW-2191 or ToE closure.

## Hard limits

- No legacy physical-role claim is transferred onto K_strict_gate.
- No alpha_geo electroweak role theorem is exported.
- No beta_tors electromagnetic/gravity hierarchy role theorem is exported.
- No beta_tors -> chi_11 theorem is exported.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
