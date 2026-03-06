# RAPORT QW-1884: NODE-STATE PARETO OOS REBALANCING

- Data UTC: 2026-03-03T00:43:58.520324+00:00
- Verdict: **PARETO_OOS_REBALANCING_PARTIAL_TRADEOFF**
- selection_mode: `feasible_min_objective`
- selected lambda_c/lambda_b: `0.2` / `0.35`

## Val Feasibility Rule
- nonboundary_rate >= 0.50
- canon_median >= 0.80
- rmse_ratio_vs_1880baseline <= 1.35

## Test Summary
- rmse median: 0.1606
- canon median: 0.9337
- nonboundary rate: 0.667

## Delta vs QW-1880
- rmse gain: -5.7721e-02
- canon gain: 4.7933e-01
- nonboundary gain: 6.6667e-01
- rmse ratio: 1.5609

## Delta vs QW-1883
- rmse gain: 3.0581e-04
- canon gain: -1.1533e-03
- nonboundary gain: 0.0000e+00

## Artifacts
- JSON: `report_qw1884_node_state_pareto_oos_rebalancing.json`
