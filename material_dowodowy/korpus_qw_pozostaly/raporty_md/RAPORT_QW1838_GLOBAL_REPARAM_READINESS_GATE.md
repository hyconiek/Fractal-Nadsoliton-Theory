# RAPORT QW-1838: GLOBAL REPARAM READINESS GATE

- Data UTC: 2026-03-02T22:30:54.042082+00:00
- Global score (QW-1835 -> now): 0.733 -> 0.900
- Delta score: 0.167
- Hard gate: **PASS**
- Readiness: **GLOBAL_CONDITIONAL_READY_UNDER_REPARAM_CRITERIA**
- Recommendation: **START_PRE_REGISTERED_JOINT_CONFIRMATORY_CAMPAIGN_PTA_GW**
- Note: Readiness is conditional on reparameterized GW criterion replacing legacy near-target criterion.

## Checks
- PTA branch readiness (1824): PASS | score=1.000 | note=SEQUENCE_BRANCH_CONDITIONAL_READY_UNDER_QUANTILE_GATING
- GW control-calibrated objective (1836): PASS | score=1.000 | note=GW_CONTROL_CALIBRATED_OBJECTIVE_SUPPORTED
- GW criterion-transition justification (1837): PASS | score=1.000 | note=GW_READY_UNDER_REPARAM_CRITERION
- Legacy transfer status (1825) [for traceability]: PARTIAL | score=0.600 | note=QUANTILE_PROTOCOL_PTA_READY_GW_BLOCKED

## Artifacts
- JSON: `report_qw1838_global_reparam_readiness_gate.json`
