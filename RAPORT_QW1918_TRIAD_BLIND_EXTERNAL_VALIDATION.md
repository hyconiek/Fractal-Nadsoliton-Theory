# RAPORT QW-1918: TRIAD BLIND EXTERNAL VALIDATION

- Data UTC: 2026-03-03T03:54:39.086306+00:00
- Triad: omega=0.088000, phi=0.890000, beta=2.000000
- Verdict: **TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG**

## Primary Dataset (blind holdout)
- path: `/home/krzysiek/Pobrane/TOE/edison/external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg/pta_v2_pairs.csv`
- pearson_corr: 0.6577 (q95 null: 0.0517, p=0.0002)
- spearman_corr: 0.3201
- rmse_gain_ratio: 0.2467 (q95 null: -0.1712, p=0.0002)
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
- JSON: `report_qw1918_triad_blind_external_validation.json`
