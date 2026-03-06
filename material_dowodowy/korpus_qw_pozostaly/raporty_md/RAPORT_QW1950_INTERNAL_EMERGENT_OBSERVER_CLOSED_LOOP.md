# RAPORT QW-1950: INTERNAL EMERGENT OBSERVER CLOSED LOOP

- Data UTC: 2026-03-03T07:07:03.627585+00:00
- Verdict: **INTERNAL_OBSERVER_CHANNEL_FAIL**

## Closed Loop Metrics
- accuracy: 0.6111
- mean reconstruction corr: 0.9349
- info bits: 4.5460
- observer loop gain: 0.0823
- observer state stability: 0.0039

## Relative Gains
- accuracy gain (closed-open): 0.0139
- info gain (closed-control): -0.0222

## Flags
- closed_accuracy_ge_min: False
- closed_info_bits_ge_min: True
- closed_mean_corr_ge_min: True
- closed_acc_gain_vs_open_ge_min: False
- closed_info_gain_vs_control_ge_min: False
- observer_loop_gain_ge_min: True
- observer_state_stability_le_max: True
- closed_loop_pass: False

## Link to Mass
- mass_strict_pass_from_qw1947: False

## Required Next Step
- REWORK_INTERNAL_OBSERVER_FEEDBACK_OPERATOR_AT_KERNEL_LEVEL

## Artifacts
- JSON: `report_qw1950_internal_emergent_observer_closed_loop.json`
