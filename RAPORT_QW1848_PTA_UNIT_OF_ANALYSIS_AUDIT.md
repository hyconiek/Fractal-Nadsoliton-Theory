# RAPORT QW-1848: PTA UNIT OF ANALYSIS AUDIT

- Data UTC: 2026-03-02T23:05:11.779944+00:00
- Cohort pairs / observed pairs: 90 / 89
- Replications: 80
- Pair obs per index (min/med/max): 10 / 22.0 / 37

## Split vs Pair
- Split-level P(rep_mean_gain>0): 0.988
- Pair-level P(pair_mean_gain>0): 0.820
- Compression gap: 0.167

## Pair-Level Inference
- k/n positives: 73/89
- one-sided lower95 for prob: 0.740
- p-value H0: prob<=0.90 -> 0.993306
- p-value H0: prob<=0.80 -> 0.374316
- Verdict: **PTA_UNIT_MISMATCH_REQUIRES_CRITERION_REDESIGN**

## Pass Flags
- split_level_prob_ge_0p9: True
- pair_level_prob_ge_0p9_inferential: False
- pair_level_prob_ge_0p8_inferential: False

## Artifacts
- JSON: `report_qw1848_pta_unit_of_analysis_audit.json`
