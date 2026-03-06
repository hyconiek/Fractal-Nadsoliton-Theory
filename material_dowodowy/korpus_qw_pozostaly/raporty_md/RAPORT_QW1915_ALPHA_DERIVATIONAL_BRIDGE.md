# RAPORT QW-1915: ALPHA DERIVATIONAL BRIDGE

- Data UTC: 2026-03-03T03:35:49.261454+00:00
- Verdict: **ALPHA_DERIVATIONAL_BRIDGE_COMPATIBLE**

## Mapping
- formula: `alpha = lambda_b / scale`
- scale: 0.05
- alpha_weighted_inv_objective: 5.690
- alpha_selected_q1891: 7.000
- alpha_grid_range_from_q1891: [2.000, 10.000]

## Empirical Reference (QW-1913)
- alpha_selected_multisplit_median: 6.000
- alpha_selected_multisplit_values: [6.0, 6.0, 6.0]
- alpha_selected_multisplit_std: 0.000

## Compatibility
- abs_diff_weighted_vs_empirical: 0.310
- abs_diff_selected_vs_empirical: 1.000
- flags: {'empirical_alpha_inside_derivational_grid': True, 'weighted_bridge_abs_diff_le_1': True, 'selected_bridge_abs_diff_le_1p5': True, 'empirical_alpha_stable_multisplit': True}

## Artifacts
- JSON: `report_qw1915_alpha_derivational_bridge.json`
