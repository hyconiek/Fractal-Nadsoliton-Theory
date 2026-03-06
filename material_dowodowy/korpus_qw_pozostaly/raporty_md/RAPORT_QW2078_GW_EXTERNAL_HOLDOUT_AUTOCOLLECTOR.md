# RAPORT QW-2078: GW EXTERNAL HOLDOUT AUTOCOLLECTOR

- Date UTC: 2026-03-04T04:03:51.630454+00:00
- Verdict: **GW_EXTERNAL_HOLDOUT_AUTOCOLLECTED**
- Input: `/home/krzysiek/Pobrane/TOE/edison/gw1831_window_features.csv`
- SHA256: `20b63aa8188e0a5a3e41fca687322ba0f66a4207f2fc307d73af95d81e3c92ce`
- Rows: `759`

## Metrics
- auc_h1l1_vs_ctrl: `0.814963521`
- adv_shared_minus_ctrl_q90: `0.310276680`
- sep_median_h1l1_minus_ctrl: `0.002055894`
- control_median_gap: `0.001288983`

## Threshold Checks
- auc_ge_min: `True`
- adv_ge_min: `True`
- sep_ge_min: `True`
- gap_le_max: `True`
- all_thresholds_pass: `True`

## Output
- observation JSON for QW-2077: `/home/krzysiek/Pobrane/TOE/edison/empirical_observations_input_qw2077.gw_autocollected.json`
- report JSON: `report_qw2078_gw_external_holdout_autocollector.json`
