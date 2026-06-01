# P2401 S1351: role-successor unlock-order certificate

Status: `PASS_ROLE_SUCCESSOR_UNLOCK_ORDER_ENUMERATED_NO_PREFIX_TRANSFER`

## Result

P2401/S1351 enumerates all six proof orders for the three role-successor atoms and measures when role claims unlock.

## Unlock order

- Best total unlock-step sum: `9`.
- Best orders: `[['alpha_geo_electroweak_role_theorem', 'beta_tors_strict_role_theorem', 'beta_power_hierarchy_successor_theorem']]`.
- Expected unlock steps: `{'legacy_weinberg_sin2_theta_role_transfer': '2', 'legacy_alpha_em_inverse_role_transfer': '8/3', 'legacy_beta_power_gravity_hierarchy_successor': '8/3', 'all_three_physical_role_transfers': '3'}`.
- First-atom summary: `{'alpha_geo_electroweak_role_theorem': {'order_count': 2, 'min_total_unlock_step_sum': 9, 'mean_total_unlock_step_sum': '19/2', 'step1_unlocked_claims_union': ['legacy_weinberg_sin2_theta_role_transfer']}, 'beta_tors_strict_role_theorem': {'order_count': 2, 'min_total_unlock_step_sum': 10, 'mean_total_unlock_step_sum': '21/2', 'step1_unlocked_claims_union': []}, 'beta_power_hierarchy_successor_theorem': {'order_count': 2, 'min_total_unlock_step_sum': 11, 'mean_total_unlock_step_sum': '11', 'step1_unlocked_claims_union': []}}`.

## Hard limits

- No role-successor atom is exported by an ordering certificate.
- No one-atom or two-atom prefix licenses role transfer or ToE closure.
- No L_total or SM/GR role-bearing term is promoted by the best order.

## Gatekeepers

`{'rg_nonduplication_audit_ran': True, 'six_orders_enumerated': True, 'unique_best_order_is_alpha_then_beta_tors_then_beta_power': True, 'best_total_unlock_step_sum_is_nine': True, 'alpha_first_has_step1_weinberg_only': True, 'beta_tors_first_unlocks_no_step1_claim': True, 'p2400_full_mask_closure_inherited': True, 'fingerprint_stable': True}`
