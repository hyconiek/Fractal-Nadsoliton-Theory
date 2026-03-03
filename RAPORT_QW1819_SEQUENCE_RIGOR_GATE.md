# RAPORT QW-1819: SEQUENCE RIGOR GATE

- Data UTC: 2026-03-02T21:49:10.290275+00:00
- Global score: 0.517
- Hard gate: **FAIL**
- Readiness: **SEQUENCE_BRANCH_NOT_READY**
- Recommendation: **PARK_AND_REDESIGN_MODELING_ASSUMPTIONS**

## Checks
- In-sample evidence (1815): PASS | score=1.000 | note=MULTISCALE_SEQUENCE_EMBEDDING_PARTIAL
- OOS predictive support (1817): PASS | score=0.931 | note=SEQUENCE_OOS_VALIDATION_PARTIAL
- Stability and dispersion control (1815/1817/1818): FAIL | score=0.139 | note=dispersion remains high
- Regime-aware stabilization (1816): FAIL | score=0.000 | note=REGIME_AWARE_SEQUENCE_MODEL_WEAK

## Inconsistencies
- In-sample dispersion criterion fails for embedding branch (std delta > 0.30).
- OOS superiority vs flat is not universal (P(M2E>flat) < 0.95).
- Robust preprocessing did not reduce OOS dispersion; variability source remains unresolved.
- Regime-aware extension degraded evidence vs base embedding model.

## Artifacts
- JSON: `report_qw1819_sequence_rigor_gate.json`
