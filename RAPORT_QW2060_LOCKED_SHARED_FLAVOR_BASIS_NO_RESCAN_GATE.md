# RAPORT QW-2060: LOCKED SHARED FLAVOR BASIS NO-RESCAN GATE

- Data UTC: 2026-03-04T02:06:49.897344+00:00
- Verdict: **LOCKED_SHARED_FLAVOR_BASIS_NO_RESCAN_PARTIAL**
- pass_count: 10/12

## Locked Basis
- p_amp/r_dist: 0.700/-0.200
- q_nu: [2, 1, 0]

## Metrics
- mass mean/max/tau-charm rel%: 12.051/34.013/14.013
- flavor CKM/PMNS mean rel%: 11.867/9.386
- GW auc/adv/sep/gap: 0.8427/0.4012/0.003868/0.003124

## Flags
- mass_mean_rel_pct_le_max: True
- mass_max_rel_pct_le_max: True
- mass_tau_charm_ratio_err_le_max: True
- ckm_mean_rel_pct_le_max: True
- pmns_mean_rel_pct_le_max: True
- gw_sep_ge_min: True
- gw_adv_ge_min: True
- gw_auc_ge_min: True
- gw_control_gap_le_max: False
- run_has_no_rescan: True
- run_has_no_refit: True
- strict_first_principles_from_kernel_only: False

## Required Next Step
- PROMOTE_TO_STRICT_FIRST_PRINCIPLES_BY_DERIVING_LOCKED_COEFFICIENTS_FROM_KERNEL_INVARIANTS

## Artifacts
- JSON: `report_qw2060_locked_shared_flavor_basis_no_rescan_gate.json`
