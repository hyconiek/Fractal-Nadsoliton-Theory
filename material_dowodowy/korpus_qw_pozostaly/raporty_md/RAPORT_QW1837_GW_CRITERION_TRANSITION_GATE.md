# RAPORT QW-1837: GW CRITERION TRANSITION GATE

- Data UTC: 2026-03-02T22:30:16.530430+00:00
- Global score: 0.950
- Hard gate: **PASS**
- Readiness: **GW_READY_UNDER_REPARAM_CRITERION**
- Recommendation: **USE_CONTROL_CALIBRATED_OBJECTIVE_AS_PRIMARY_GW_GATE**

## Checks
- Legacy near-target criterion invalidity (1829): PASS | score=1.000 | note=GW_NEAR_TARGET_REQUIRES_STRUCTURAL_REDESIGN
- New control-calibrated criterion support (1836): PASS | score=1.000 | note=GW_CONTROL_CALIBRATED_OBJECTIVE_SUPPORTED
- Objective improvement magnitude (1836 vs raw): PASS | score=1.000 | note=delta_auc=0.098, delta_adv=0.368, gap_red=0.954
- Continuity with multi-detector consistency track: PASS | score=0.800 | note=GW_MULTI_DETECTOR_PARTIAL_CONTROL_MISMATCH

## Artifacts
- JSON: `report_qw1837_gw_criterion_transition_gate.json`
