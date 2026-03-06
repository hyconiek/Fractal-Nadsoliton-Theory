# RAPORT QW-1880: NODE-STATE STRICT OOS

- Data UTC: 2026-03-03T00:35:18.218241+00:00
- Verdict: **NODE_STATE_STRICT_OOS_FAIL**
- split train/val/test: 16/6/6
- locked nodes model: M_A_2_5_8_11_or_2plus3n -> [2, 5, 8, 11]
- best lambdas: lambda_c=0.2, lambda_p=0.0

## Summaries
- train rmse/canon/nonboundary: 0.0647 / 0.0000 / 1.000
- val rmse/canon/nonboundary: 0.1217 / 0.9019 / 0.000
- test rmse/canon/nonboundary: 0.1029 / 0.4543 / 0.000
- test baseline rmse/canon/nonboundary: 0.1716 / 0.0000 / 0.000

## Test Delta vs Baseline
- rmse median gain: 6.8739e-02
- canon median gain: 4.5433e-01
- nonboundary gain: 0.0000e+00

## Artifacts
- JSON: `report_qw1880_node_state_strict_oos.json`
