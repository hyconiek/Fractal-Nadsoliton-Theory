# RAPORT QW-1909: EXTERNAL REBUILD COMPARISON

- Data UTC: 2026-03-03T02:36:36.268877+00:00
- Verdict: **BEST_CANDIDATE_PARTIAL_OR_EXTERNALITY_BLOCKED**
- readiness: **EMPIRICAL_CLOSURE_NOT_READY**

## Candidates
- `/home/krzysiek/Pobrane/TOE/edison/external_confirmatory_v2/confirmatory_dataset_external_source_rebuild`: hard=FAIL, pta_all=False, gw_all=False, externality_ok=False, pta_prob=0.467, gw_auc=0.724
- `/home/krzysiek/Pobrane/TOE/edison/external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg`: hard=PARTIAL, pta_all=False, gw_all=True, externality_ok=False, pta_prob=0.467, gw_auc=0.942

## Best Candidate
- dir: `/home/krzysiek/Pobrane/TOE/edison/external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg`
- hard_gate: PARTIAL
- externality_ok: False
- PTA mean gain: -0.000020
- PTA prob positive: 0.467
- GW mean AUC: 0.942
- GW mean adv: 0.757

## Artifacts
- JSON: `report_qw1909_external_rebuild_comparison.json`
