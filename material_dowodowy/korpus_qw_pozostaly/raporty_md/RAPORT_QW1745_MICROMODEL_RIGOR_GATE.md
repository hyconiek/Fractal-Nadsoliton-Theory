# RAPORT QW-1745: MICROMODEL RIGOR GATE

- Data UTC: 2026-03-02T17:16:56.497640+00:00
- Global score: 0.400
- Hard gate: **FAIL**
- Readiness: **MICROMODEL_ITERATION_OPEN**

## Checks
- Signed dynamic derivation (1739): PASS | score=0.800 | note=SIGNED_MICROMODEL_DERIVATION_SUPPORTED_BUT_SHIFTED
- Identifiability audit (1740): FAIL | score=0.300 | note=IDENTIFIABILITY_WEAK
- Constrained global fit (1741): PASS | score=0.750 | note=CONSTRAINED_GLOBAL_DERIVATION_NOT_CLOSED
- Profile-likelihood identifiability (1742): FAIL | score=0.000 | note=PROFILE_IDENTIFIABILITY_WEAK
- Oscillatory cohort derivation (1743): FAIL | score=0.550 | note=OSCILLATORY_COHORT_DERIVATION_NOT_CLOSED
- Oscillatory cohort identifiability (1744): FAIL | score=0.000 | note=OSCILLATORY_IDENTIFIABILITY_WEAK

## Artifacts
- JSON: `report_qw1745_micromodel_rigor_gate.json`
