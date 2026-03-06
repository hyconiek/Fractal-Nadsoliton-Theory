# RAPORT QW-1830: GLOBAL EMPIRICAL STATUS GATE

- Data UTC: 2026-03-02T22:10:22.026202+00:00
- Global score: 0.444
- Hard gate: **PARTIAL**
- Readiness: **GLOBAL_PARTIAL_PTA_ONLY**
- Recommendation: **RUN_PTA_TRACK_CONTINUATION_AND_PARALLEL_GW_REDESIGN_PROGRAM**

## Checks
- PTA quantile-gated branch: PASS | score=1.000 | note=SEQUENCE_BRANCH_CONDITIONAL_READY_UNDER_QUANTILE_GATING
- Cross-domain transfer status: PARTIAL | score=0.600 | note=QUANTILE_PROTOCOL_PTA_READY_GW_BLOCKED
- GW branch readiness: FAIL | score=0.176 | note=GW_BRANCH_REDESIGN_REQUIRED
- GW near-target feasibility: FAIL | score=0.000 | note=GW_NEAR_TARGET_REQUIRES_STRUCTURAL_REDESIGN

## Next Steps
- QW-1831: event-windowed GW coherent feature extractor (chirp-conditioned).
- QW-1832: GW shared-vs-control objective with explicit near-target prevalence optimization.
- QW-1833: GW multi-detector consistency gate on H1-L1, H1-V1, L1-V1 with quantile scoring.

## Artifacts
- JSON: `report_qw1830_global_empirical_status_gate.json`
