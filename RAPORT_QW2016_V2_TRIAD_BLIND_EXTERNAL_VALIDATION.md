# RAPORT QW-2016: TRIAD BLIND EXTERNAL VALIDATION

- Data UTC: 2026-03-04T07:35:04.040096+00:00
- Triad: omega=0.088000, phi=0.890000, beta=2.000000
- Verdict: **TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG**

## Primary Dataset (blind holdout)
- path: `/home/krzysiek/Pobrane/TOE/edison/external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv`
- pearson_corr: 0.0531 (q95 null: 0.0369, p=0.009198)
- spearman_corr: 0.0308
- rmse_gain_ratio: 0.0003 (q95 null: -0.0013, p=0.009198)
- all_pass: True

## Stress Dataset (blind holdout)
- path: `/home/krzysiek/Pobrane/TOE/edison/external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv`
- pearson_corr: 0.2932 (q95 null: 0.0519, p=0.0002)
- spearman_corr: 0.1696
- rmse_gain_ratio: 0.0437 (q95 null: -0.0227, p=0.0002)
- all_pass: True

## Gate Flags
- primary_all_pass: True
- stress_soft_pass: True

## Artifacts
- JSON: `report_qw2016_v2_triad_blind_external_validation.json`
