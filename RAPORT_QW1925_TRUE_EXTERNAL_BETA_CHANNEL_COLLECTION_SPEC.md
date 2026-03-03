# RAPORT QW-1925: TRUE EXTERNAL BETA-CHANNEL COLLECTION SPEC

- Data UTC: 2026-03-03T19:06:53.947886+00:00
- Verdict: **TRUE_EXTERNAL_BETA_CHANNEL_COLLECTION_SPEC_READY**
- selected_beta_observable: `B7_local_resid_std`

## Minimal Signal Targets
- holdout_effect_beta >= 0.6
- holdout_effect_omega <= 0.25
- holdout_contrast >= 0.35
- holdout_contrast_boot_q05 >= 0.2

## Power Targets
- conservative contrast target: 0.7904
- conservative sigma: 0.1200
- n_holdout min (80% power): 400
- n_holdout min (90% power): 500
- n_total_pairs recommended: 1200

## Required Files
- manifest_beta_channel.json
- beta_channel_pairs.csv
- intervention_events.csv
- protocol_freeze.json

## Artifacts
- JSON: `report_qw1925_true_external_beta_channel_collection_spec.json`
