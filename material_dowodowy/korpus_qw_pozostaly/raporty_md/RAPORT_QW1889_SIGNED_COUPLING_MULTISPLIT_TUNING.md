# RAPORT QW-1889: SIGNED-COUPLING MULTISPLIT TUNING

- Data UTC: 2026-03-03T00:55:20.034854+00:00
- Verdict: **SIGNED_COUPLING_TUNED_HOLDOUT_WEAK**
- selection_mode: `feasible_max_objective`
- selected lambda_c/lambda_b: `0.2` / `0.35`

## Holdout Summary (signed - unsigned)
- success_rate: 0.462
- rmse_gain_median: 0.0162
- canon_gain_median: -0.0660
- nonboundary_gain_median: 0.3333
- rmse_gain_q25: 0.0101
- canon_gain_q25: -0.4453
- nonboundary_gain_q25: 0.1667

## Success Rule
- rmse_gain >= 0.005
- canon_gain >= -0.10
- nonboundary_gain >= 0.0

## Artifacts
- JSON: `report_qw1889_signed_coupling_multisplit_tuning.json`
