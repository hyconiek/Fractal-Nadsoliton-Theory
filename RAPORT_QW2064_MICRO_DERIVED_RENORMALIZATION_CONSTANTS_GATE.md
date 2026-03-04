# RAPORT QW-2064: MICRO-DERIVED RENORMALIZATION CONSTANTS GATE

- Data UTC: 2026-03-04T02:23:16.916709+00:00
- Verdict: **MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS_WITH_WIDE_CI_WARNING**
- pass_count: 8/8
- ci_warning: True

## Targets (from frozen kernel)
- z_beta_target: 100.000000
- delta_eta_target: 0.800000

## Micro-derived (QW-2048)
- z_beta_median: 114.739580
- delta_eta_median: 1.125000
- n_rows_total / n_bins: 342 / 17
- phase_min_median: 0.852165
- joint_target_in_ci95_fraction: 0.823529

## Deviations
- abs_log(z_beta_micro/z_beta_target): 0.137495
- abs(delta_eta_micro-delta_eta_target): 0.325000

## Dispersion Diagnostics
- z_beta log-IQR: 3.132919
- delta_eta IQR: 0.350000

## Flags
- qw2048_pointwise_pass: True
- n_rows_total_ge_min: True
- n_bins_ge_min: True
- phase_min_median_ge_min: True
- joint_target_in_ci95_fraction_ge_min: True
- z_beta_logratio_to_target_le_max: True
- delta_eta_abs_err_le_max: True
- no_sector_retune_between_micro_and_kernel: True

## Artifacts
- JSON: `report_qw2064_micro_derived_renormalization_constants_gate.json`
