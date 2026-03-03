# RAPORT QW-1839: JOINT CONFIRMATORY PREREG PROTOCOL

- Data UTC: 2026-03-02T22:32:27.188597+00:00
- Verdict: **JOINT_CONFIRMATORY_PREREG_FROZEN**
- Protocol SHA256: `b9cf21d3d32508e95c6f7cef2e8a953a12f6c7ea7732b6288605c518ab5db5af`

## PTA Thresholds
- mean_quantile_gain >= 0.04
- prob_quantile_gain_positive >= 0.9
- std_quantile_gain <= 0.035

## GW Thresholds
- calibrated_mean_auc >= 0.9
- calibrated_mean_adv >= 0.6
- calibrated_mean_control_gap <= 0.0005
- calibrated_prob_adv_positive >= 0.9

## Anti-Leakage
- Offsets/calibrations train-only.
- Split rule frozen (window_idx mod 5).
- No post-hoc threshold tuning.

## Artifacts
- JSON: `report_qw1839_joint_confirmatory_prereg_protocol.json`
