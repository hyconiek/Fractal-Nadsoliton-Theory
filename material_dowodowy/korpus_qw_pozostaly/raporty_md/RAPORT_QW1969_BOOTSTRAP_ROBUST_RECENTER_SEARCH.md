# RAPORT QW-1969: BOOTSTRAP ROBUST RECENTER SEARCH

- Data UTC: 2026-03-03T08:34:12.334271+00:00
- Verdict: **INSUFFICIENT_BOOTSTRAP_ROBUSTNESS**

## Baseline vs Best Recentered
- baseline triad bootstrap pass (5000 from QW-1968): 71.92%
- best recentered triad bootstrap pass (5000): 71.94%
- deterministic candidate pool (triad pass): 55544 / 60001

## Best Recentered Deterministic Metrics
- CKM/PMNS mean rel%: 10.242/11.436
- GW auc/adv/sep/gap: 0.8238/0.3221/0.002721/0.002308

## Local Deterministic Pass Around Best Recentered
- radius=0.010: 100.00% (8000/8000)
- radius=0.020: 100.00% (8000/8000)
- radius=0.050: 86.14% (6891/8000)

## Required Next Step
- ADD_STRUCTURAL_GW_CONTROL_TERM_WITH_SHARED_ORIGIN_AND_REPEAT

## Artifacts
- JSON: `report_qw1969_bootstrap_robust_recenter_search.json`
