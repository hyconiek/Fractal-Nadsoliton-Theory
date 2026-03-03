# RAPORT QW-1933: REPARAM HASH-SPLIT ROBUSTNESS

- Data UTC: 2026-03-03T05:39:50.230764+00:00
- Triad (eta branch): omega=0.373414, phi=-1.310234, beta=0.615938, eta=2.80
- Verdict: **REPARAM_HASH_SPLIT_ROBUSTNESS_PASS**

## Primary Summary
- pass_rate: 0.8750 (21/24)
- corr median [q10, q90]: 0.0706 [0.0495, 0.0984]
- gain median [q10, q90]: 0.0024 [0.0002, 0.0035]

## Stress Summary
- pass_rate: 1.0000 (24/24)
- corr median [q10, q90]: 0.3177 [0.3007, 0.3514]
- gain median [q10, q90]: 0.0508 [0.0446, 0.0627]

## Strict Flags
- primary_pass_rate_ge_0p70: True
- stress_pass_rate_ge_0p90: True
- primary_corr_median_positive: True
- primary_gain_median_positive: True
- stress_corr_median_positive: True
- stress_gain_median_positive: True

## Required Next Step
- INTEGRATE_QW1932_QW1933_IN_STRICT_CLOSURE_GATE

## Artifacts
- JSON: `report_qw1933_reparam_hash_split_robustness.json`
