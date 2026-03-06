# RAPORT QW-1825: QUANTILE PROTOCOL TRANSFER GATE

- Data UTC: 2026-03-02T22:04:47.518309+00:00
- Global score: 0.741
- Hard gate: **PARTIAL**
- Readiness: **QUANTILE_PROTOCOL_PTA_READY_GW_BLOCKED**
- Recommendation: **RUN_GW_SPECIFIC_REDESIGN_BEFORE_JOINT_CAMPAIGN**
- PTA ready: **True**
- GW ready: **False**

## Checks
- PTA quantile predictive branch (1823): PASS | score=1.000 | note=QUANTILE_SCORE_OOS_SUPPORTED
- Quantile-gated readiness decision (1824): PASS | score=1.000 | note=SEQUENCE_BRANCH_CONDITIONAL_READY_UNDER_QUANTILE_GATING
- PTA operational campaign robustness (1786): PASS | score=0.966 | note=EMPIRICAL_CAMPAIGN_ROBUSTNESS_CONFIRMED
- GW transfer feasibility (1725+1726): FAIL | score=0.000 | note=GW_CROSS_HURST_ANOMALY_NOT_ROBUST / FIN_023_TO_031_PROJECTION_NOT_SUPPORTED

## Blockers
- GW branch fails strict robustness and projection criteria (QW-1725/QW-1726).

## Artifacts
- JSON: `report_qw1825_quantile_protocol_transfer_gate.json`
