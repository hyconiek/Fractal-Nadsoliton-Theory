# RAPORT QW-1812: DYNAMIC PROXY BRANCH GATE

- Data UTC: 2026-03-02T21:18:02.778371+00:00
- Global score: 0.200
- Hard gate: **FAIL**
- Readiness: **DYNAMIC_PROXY_BRANCH_NOT_READY**
- Recommendation: **MOVE_TO_SEQUENCE_LEVEL_DYNAMIC_MODELS**

## Checks
- Drift proxy model (1808): FAIL | score=0.200 | note=DYNAMIC_DRIFT_REGIME_WEAK
- Feature scan proxies (1809): FAIL | score=0.000 | note=DYNAMIC_FEATURE_SCAN_WEAK
- Triad dynamic proxy model (1811): FAIL | score=0.400 | note=DYNAMIC_TRIAD_MODEL_WEAK

## Artifacts
- JSON: `report_qw1812_dynamic_proxy_branch_gate.json`
