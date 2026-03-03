# RAPORT QW-1906: EXTERNAL DATA COLLECTION SPEC

- Data UTC: 2026-03-03T02:04:00.646496+00:00
- Verdict: **EXTERNAL_DATA_COLLECTION_SPEC_READY**

## Minimal External Targets (PTA)
- n_pairs recommended min: 1200
- n_pairs preferred: 2000
- mean_pair_mean_gain >= 0.04
- prob_pair_mean_gain_positive >= 0.667
- lower95_prob_positive >= 0.6

## Feature-Signal Proxy Targets
- alpha_required_for_80pct_power: 6.0
- hxy-feature slope baseline -> target: 0.0079 -> 0.2744
- hxy-feature corr baseline -> target: 0.1751 -> 0.9724
- delta_hxy |mean| needed: 0.2162
- delta_hxy positive fraction: 0.428

## Artifacts
- JSON: `report_qw1906_external_data_collection_spec.json`
