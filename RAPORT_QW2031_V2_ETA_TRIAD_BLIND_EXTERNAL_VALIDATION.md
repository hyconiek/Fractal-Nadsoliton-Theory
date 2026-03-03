# RAPORT QW-2031: V2 ETA TRIAD BLIND EXTERNAL VALIDATION

- Data UTC: 2026-03-03T19:48:59.411505+00:00
- Kernel: omega=0.233750, phi=-0.137500, beta=1.170000, eta=1.800000
- Readiness: **READY_FOR_COMBINED_CONFIRMATORY_GATE**
- Verdict: **ETA_TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG**

## Primary Dataset (blind holdout)
- path: `/home/krzysiek/Pobrane/TOE/edison/external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv`
- pearson_corr: 0.0546 (q95 null: 0.0367, p=0.007598)
- spearman_corr: 0.0308
- rmse_gain_ratio: 0.0006 (q95 null: -0.0012, p=0.007598)
- all_pass: True

## Stress Dataset (blind holdout)
- path: `/home/krzysiek/Pobrane/TOE/edison/external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv`
- pearson_corr: 0.3300 (q95 null: 0.0512, p=0.0002)
- spearman_corr: 0.1696
- rmse_gain_ratio: 0.0558 (q95 null: -0.0318, p=0.0002)
- all_pass: True

## Gate Flags
- primary_all_pass: True
- stress_soft_pass: True

## Required Next Step
- RUN_QW2032_COMBINED_CONFIRMATORY_GATE

## Artifacts
- JSON: `report_qw2031_v2_eta_triad_blind_external_validation.json`
