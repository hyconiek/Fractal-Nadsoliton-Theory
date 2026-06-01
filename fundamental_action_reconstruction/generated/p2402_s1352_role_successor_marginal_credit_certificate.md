# P2402 S1352: role-successor marginal-credit certificate

Status: `PASS_ROLE_SUCCESSOR_MARGINAL_CREDIT_PRIORITIZES_ALPHA_NO_PREFIX_LICENSE`

## Result

P2402/S1352 computes role-local marginal credit over all six proof orders for the three role-successor atoms.

## Marginal credit

- Ranking over all claims: `['alpha_geo_electroweak_role_theorem', 'beta_tors_strict_role_theorem', 'beta_power_hierarchy_successor_theorem']`.
- Ranking over physical claims: `['alpha_geo_electroweak_role_theorem', 'beta_tors_strict_role_theorem', 'beta_power_hierarchy_successor_theorem']`.
- Atom credit: `{'alpha_geo_electroweak_role_theorem': {'mean_marginal_claim_count': '11/6', 'mean_physical_marginal_claim_count': '3/2', 'marginal_count_samples_by_order': [1, 1, 2, 3, 1, 3], 'physical_marginal_count_samples_by_order': [1, 1, 2, 2, 1, 2], 'claim_credit_breakdown': {'legacy_weinberg_sin2_theta_role_transfer': '1', 'legacy_alpha_em_inverse_role_transfer': '1/2', 'legacy_beta_power_gravity_hierarchy_successor': '0', 'all_three_physical_role_transfers': '1/3'}}, 'beta_tors_strict_role_theorem': {'mean_marginal_claim_count': '4/3', 'mean_physical_marginal_claim_count': '1', 'marginal_count_samples_by_order': [1, 3, 0, 0, 3, 1], 'physical_marginal_count_samples_by_order': [1, 2, 0, 0, 2, 1], 'claim_credit_breakdown': {'legacy_weinberg_sin2_theta_role_transfer': '0', 'legacy_alpha_em_inverse_role_transfer': '1/2', 'legacy_beta_power_gravity_hierarchy_successor': '1/2', 'all_three_physical_role_transfers': '1/3'}}, 'beta_power_hierarchy_successor_theorem': {'mean_marginal_claim_count': '5/6', 'mean_physical_marginal_claim_count': '1/2', 'marginal_count_samples_by_order': [2, 0, 2, 1, 0, 0], 'physical_marginal_count_samples_by_order': [1, 0, 1, 1, 0, 0], 'claim_credit_breakdown': {'legacy_weinberg_sin2_theta_role_transfer': '0', 'legacy_alpha_em_inverse_role_transfer': '0', 'legacy_beta_power_gravity_hierarchy_successor': '1/2', 'all_three_physical_role_transfers': '1/3'}}}`.

## Hard limits

- No marginal-credit score exports a role-successor atom.
- No marginal-credit score licenses a one-atom or two-atom prefix for role transfer.
- No L_total or SM/GR role-bearing term is promoted by this prioritization metric.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'six_orders_enumerated': True, 'alpha_total_credit_is_11_over_6': True, 'beta_tors_total_credit_is_4_over_3': True, 'beta_power_total_credit_is_5_over_6': True, 'alpha_is_top_total_and_physical_credit': True, 'p2401_best_order_inherited': True, 'p2400_full_role_mask_still_required': True, 'fingerprint_stable': True}`
