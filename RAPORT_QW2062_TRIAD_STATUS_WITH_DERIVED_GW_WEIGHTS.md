# RAPORT QW-2062: TRIAD STATUS WITH DERIVED GW WEIGHTS

- Data UTC: 2026-03-04T02:09:16.856878+00:00
- Verdict: **TRIAD_PHYSICAL_THRESHOLDS_PASS_WITH_LOCKED_FLAVOR_AND_DERIVED_GW_WEIGHTS**
- pass_count: 10/11
- physical_pass: True

## Metrics
- mass mean/max/tau-charm rel%: 12.051/34.013/14.013
- flavor CKM/PMNS mean rel%: 11.867/9.386
- GW auc/adv/sep/gap: 0.8150/0.3103/0.002056/0.001289

## Derived GW Weights
- {'w_max_abs_corr': 0.0, 'w_mean_abs_corr': 0.6424131897825203, 'w_corr_at_0ms': 0.1834822641681202, 'w_corr_at_10ms': 0.17410454604935954, 'decay_ratio': 0.10006511151920049}

## Flags
- mass_mean_rel_pct_le_max: True
- mass_max_rel_pct_le_max: True
- mass_tau_charm_ratio_err_le_max: True
- ckm_mean_rel_pct_le_max: True
- pmns_mean_rel_pct_le_max: True
- gw_sep_ge_min: True
- gw_adv_ge_min: True
- gw_auc_ge_min: True
- gw_control_gap_le_max: True
- gw_weights_derived_no_scan: True
- strict_first_principles_from_kernel_only: False

## Required Next Step
- DERIVE_LOCKED_FLAVOR_BASIS_FROM_KERNEL_INVARIANTS_TO_PROMOTE_TO_STRICT_FIRST_PRINCIPLES

## Artifacts
- JSON: `report_qw2062_triad_status_with_derived_gw_weights.json`
