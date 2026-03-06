# RAPORT QW-1962: NONCIRCULAR BRANCH UNIFIED TRIAD GATE

- Data UTC: 2026-03-03T08:06:37.787781+00:00
- Verdict: **NONCIRCULAR_UNIFIED_TRIAD_FAIL**

## Branch Input
- q_assignment: legacy_fibonacci
- gamma_used: 1.506667
- delta_info_used: 0.177027
- shared params p/r/phase: 0.9575/0.4009/1.7857

## Mass
- mean/max rel err: 12.051% / 34.013%
- tau/charm ratio pred/exp/error: 1.2031/1.3991/14.013%

## Flavor
- CKM mean rel err: 48.826%
- PMNS mean rel err: 30.062%

## GW
- auc/adv/sep/gap: 0.8426/0.3972/0.003644/0.002958

## Flags
- mass_mean_rel_pct_le_max: True
- mass_max_rel_pct_le_max: True
- mass_tau_charm_ratio_err_le_max: True
- ckm_mean_rel_pct_le_max: False
- pmns_mean_rel_pct_le_max: False
- gw_sep_ge_min: True
- gw_adv_ge_min: True
- gw_auc_ge_min: True
- gw_control_gap_le_max: False

## Required Next Step
- KEEP_NONCIRCULAR_MASS_BRANCH_AND_REWORK_SHARED_FLAVOR_GW_OPERATOR

## Artifacts
- JSON: `report_qw1962_noncircular_branch_unified_triad_gate.json`
