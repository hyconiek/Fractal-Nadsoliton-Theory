# RAPORT QW-1886: NODE-STATE PHASE-7 GATE

- Data UTC: 2026-03-03T00:47:35.797159+00:00
- Verdict: **NODE_STATE_PHASE7_GATE_COMPLETE**
- readiness: **PHASE7_PARTIAL_SINGLE_SPLIT_ONLY_REQUIRES_MODEL_REFORMULATION**
- hard_gate: **FAIL**
- global_score: 0.443

## Components
- QW-1882 strict OOS: score=0.433, pass=False
- QW-1884 pareto: score=0.694, soft_pass=True
- QW-1885 robustness: score=0.294, pass=False

## Key Metrics
- 1884 test: canon=0.934, nonboundary=0.667, rmse_ratio_vs_1880=1.561
- 1885 multisplit: success_rate=0.160, rmse_penalty_median=0.084, nonboundary_gain_median=0.333

## Next Required Step
- QW-1887_SIGNED_COUPLING_MICRODYNAMICS_REBUILD

## Artifacts
- JSON: `report_qw1886_node_state_phase7_gate.json`
