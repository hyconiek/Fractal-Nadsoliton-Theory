# RAPORT QW-1828: GW TARGET THRESHOLD MAP

- Data UTC: 2026-03-02T22:08:23.650283+00:00
- GW readiness input: GW_BRANCH_REDESIGN_REQUIRED
- Missing targets: 8
- Already met targets: 1
- Critical targets: 5
- Major targets: 4

## Targets
- prob_shared_beats_control_near_target: current=0, target=0.7, direction=up, factor_needed=inf, priority=critical, status=missing
- phase_null_p_lower: current=1, target=0.01, direction=down, factor_needed=100.000, priority=critical, status=missing
- shift_null_p_lower: current=0.801653, target=0.01, direction=down, factor_needed=80.165, priority=critical, status=missing
- corr_at_plus_10ms: current=0.0016325, target=0.02, direction=up, factor_needed=12.251, priority=critical, status=missing
- abs_diff_H1L1_vs_H1V1: current=0.446043, target=0.15, direction=down, factor_needed=2.974, priority=critical, status=missing
- mean_abs_effect_size_d: current=0.178723, target=0.35, direction=up, factor_needed=1.958, priority=major, status=missing
- mean_nonoverlap_frac: current=0.245342, target=0.35, direction=up, factor_needed=1.427, priority=major, status=missing
- length_spread: current=0.0990116, target=0.08, direction=down, factor_needed=1.238, priority=major, status=missing
- window_std: current=0.0125572, target=0.05, direction=down, factor_needed=0.251, priority=major, status=met

## Artifacts
- JSON: `report_qw1828_gw_target_threshold_map.json`
