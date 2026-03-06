# RAPORT QW-1786: CAMPAIGN ROBUSTNESS GATE

- Data UTC: 2026-03-02T19:03:07.989461+00:00
- Global score: 0.966
- Hard gate: **PASS**
- Readiness: **EMPIRICAL_CAMPAIGN_ROBUSTNESS_CONFIRMED**

## Checks
- Campaign launch status (1779): PASS | score=0.895 | note=EMPIRICAL_CAMPAIGN_STRONGLY_ON_TRACK
- Operational cohort readiness (1781): PASS | score=1.000 | note=COHORT_OPERATIONAL_READY
- Low-coverage stress robustness (1782+1784): PASS | score=0.968 | note=score_1782=0.975, score_1784=0.961
- High-coverage replication stability (1785): PASS | score=1.000 | note=HIGH_COVERAGE_REPLICATION_SUPPORTED

## Stress Details
- score_1782: 0.975
- score_1784: 0.961
- mean_stress_score: 0.968

## Artifacts
- JSON: `report_qw1786_campaign_robustness_gate.json`
