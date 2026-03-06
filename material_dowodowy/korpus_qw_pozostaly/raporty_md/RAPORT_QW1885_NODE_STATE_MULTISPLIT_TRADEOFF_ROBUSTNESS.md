# RAPORT QW-1885: NODE-STATE MULTISPLIT TRADEOFF ROBUSTNESS

- Data UTC: 2026-03-03T00:46:01.035797+00:00
- Verdict: **MULTISPLIT_TRADEOFF_NOT_ROBUST**
- Splits: 25 (seed 188000..188024)
- lambda_c/lambda_b(control->treatment): 0.2 / 0.0->0.35

## Global Summary
- success_rate: 0.160
- rmse_penalty_median: 0.0843
- rmse_penalty_q75: 0.1254
- canon_gain_median: 0.4398
- nonboundary_gain_median: 0.3333
- nonboundary_gain_q25: 0.1667

## Split Success Rule
- nonboundary_gain >= 0.30
- canon_gain >= -0.05
- rmse_penalty <= 0.07

## Artifacts
- JSON: `report_qw1885_node_state_multisplit_tradeoff_robustness.json`
