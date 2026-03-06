# RAPORT QW-2055: STRICT FIRST-PRINCIPLES TRIAD CLOSURE GATE

- Data UTC: 2026-03-04T01:45:00.215544+00:00
- Verdict: **STRICT_FIRST_PRINCIPLES_TRIAD_CLOSURE_FAIL**
- pass_count: 11/13

## Frozen Inputs
- kernel (omega/phi/beta/eta): 0.185750/0.162500/1.000000/1.800000
- q_assignment (QW-1961): legacy_fibonacci
- gamma: 1.506667
- delta_info: 0.177027
- deterministic shared params p/r/phase: 1.034539/0.824205/1.291826

## Metrics
- mass mean/max/tau-charm rel%: 12.051/34.013/14.013
- flavor CKM/PMNS mean rel%: 197.086/14.060
- GW auc/adv/sep/gap: 0.8320/0.3419/0.003108/0.002611

## Flags
- mass_mean_rel_pct_le_max: True
- mass_max_rel_pct_le_max: True
- mass_tau_charm_ratio_err_le_max: True
- ckm_mean_rel_pct_le_max: False
- pmns_mean_rel_pct_le_max: True
- gw_sep_ge_min: True
- gw_adv_ge_min: True
- gw_auc_ge_min: True
- gw_control_gap_le_max: False
- protocol_no_fit: True
- protocol_no_sector_retune: True
- protocol_kernel_source_qw2049: True
- protocol_mass_chain_source_qw1961: True

## Ranked Gaps (Descending)
- ckm_mean_rel_pct_over: 182.085761
- gw_control_gap_over: 0.000111
- mass_mean_rel_pct_over: 0.000000
- mass_max_rel_pct_over: 0.000000
- tau_charm_ratio_rel_pct_over: 0.000000
- pmns_mean_rel_pct_over: 0.000000
- gw_sep_under: 0.000000
- gw_adv_under: 0.000000
- gw_auc_under: 0.000000

## Required Next Step
- REWORK_FLAVOR_OPERATOR_FROM_KERNEL_INVARIANTS_WITHOUT_INTRODUCING_NEW_FREE_PARAMETERS

## Artifacts
- JSON: `report_qw2055_strict_first_principles_triad_closure_gate.json`
