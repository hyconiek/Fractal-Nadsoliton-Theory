# RAPORT QW-1843: PTA THRESHOLD INFERENCE RIGOR

- Data UTC: 2026-03-02T22:38:55.715093+00:00
- n replications: 14
- Observed mean/std/prob+: 0.062399 / 0.026698 / 1.000
- Thresholds mean>=/std<=/prob>=: 0.040 / 0.035 / 0.900
- Mean CI95: [0.048746, 0.075655]
- Std CI95: [0.016826, 0.033768]
- Prob-inference: p-value(H0:p<=thr)=0.228768, one-sided lower95=0.807364
- n_min(all-positive) for alpha=0.05 vs prob-threshold: 29
- Verdict: **PTA_POINT_PASS_BUT_PROBABILITY_UNDERPOWERED**

## Pass Flags
- mean_threshold_with_ci: True
- std_threshold_with_ci: True
- prob_threshold_point: True
- prob_threshold_inferential: False

## Artifacts
- JSON: `report_qw1843_pta_threshold_inference_rigor.json`
