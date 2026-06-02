# P2446 S1396: strict pointwise rank-lift selector obstruction certificate

Status: `PASS_STRICT_POINTWISE_RANK_LIFT_SELECTOR_OBSTRUCTION_NO_SOURCE_THEOREM`

## Finite facts

- Scan point count: `2001`.
- `d=1` normalized volume: `0.020925092942762662`.
- Global best point: `{'d': 0.785, 'k_strict_value': 0.5786073656846086, 'normalized_rank_lift_volume': 0.022155467429392502, 'robust_rank_lift_candidate': True}`.
- Local best point: `{'d': 0.785, 'k_strict_value': 0.5786073656846086, 'normalized_rank_lift_volume': 0.022155467429392502, 'robust_rank_lift_candidate': True}`.
- Robust pointwise count in local window: `201`.
- `d=1` is global maximum: `False`.
- `d=1` is local maximum: `False`.

## Hard limits

- Pointwise rank-lift conditioning does not select a unique evaluation coordinate d_ref.
- The stable K_at_d_1 row remains a candidate because conditioning alone neither proves strict observable/source admissibility nor supplies a lawful gauge slice.
- No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.

## Next honest step

Prove a strict point-coordinate selector, observable/source theorem, or gauge-slice theorem before using any pointwise K_strict(d_ref) row as a supplemental coefficient source.

## Gatekeepers

`{'rg_audit_ran': True, 'p2445_best_inherited': True, 'scan_grid_complete': True, 'd1_remains_robust': True, 'd1_not_global_maximum': True, 'd1_not_local_maximum': True, 'positive_global_selector_gap': True, 'positive_local_selector_gap': True, 'many_robust_pointwise_alternatives': True, 'no_pointwise_selector_export': True, 'no_observable_source_export': True, 'no_gauge_fixing_export': True, 'no_value_generator_export': True, 'no_qw2191_discharge': True, 'no_ltotal_export': True, 'no_toe_export': True, 'fingerprint_stable': True}`
