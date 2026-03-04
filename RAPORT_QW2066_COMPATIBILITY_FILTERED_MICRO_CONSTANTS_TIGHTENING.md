# RAPORT QW-2066: COMPATIBILITY-FILTERED MICRO CONSTANTS TIGHTENING

- Data UTC: 2026-03-04T02:26:05.450664+00:00
- Verdict: **COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING_PASS**
- pass_count: 6/6
- selection_mode: compatibility_constrained_min_dispersion
- tightened_warning_resolved: True

## Baseline vs Selected
- baseline z_beta_log_iqr: 3.132919
- selected z_beta_log_iqr: 2.123599
- baseline delta_eta_iqr: 0.350000
- selected delta_eta_iqr: 0.350000
- selected n_bins: 5

## Selected Medians
- z_beta_q50: 109.761471 (target 100.000000)
- delta_eta_q50: 0.956250 (target 0.800000)
- abs_log_ratio_z_beta: 0.093139
- abs_err_delta_eta: 0.156250

## Selected Filter
- q_n / q_phase / q_rmse: 0.20 / 0.40 / 0.70
- n_min / phase_min / rmse_max: 10.200 / 0.848966 / 0.100441

## Flags
- n_selected_ge_min: True
- z_beta_abs_log_ratio_to_target_le_max: True
- delta_eta_abs_err_to_target_le_max: True
- z_beta_log_iqr_le_max: True
- delta_eta_iqr_le_max: True
- phase_selected_median_ge_min: True

## Artifacts
- JSON: `report_qw2066_compatibility_filtered_micro_constants_tightening.json`
