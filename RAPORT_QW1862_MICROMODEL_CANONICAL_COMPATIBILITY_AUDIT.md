# RAPORT QW-1862: MICROMODEL CANONICAL COMPATIBILITY AUDIT

- Data UTC: 2026-03-03T00:00:00.466473+00:00
- Verdict: **MICROMODEL_CANONICAL_COMPATIBILITY_FAIL**

## Target Tuple
- omega target: 0.785398 (QW-1861 winner)
- phi target: 0.523599 (QW-1861 winner)
- beta target: 0.010000 (QW-1861 winner)

## Source Compatibility
- QW-1739: final_score=0.0000, raw_score=0.0000, fit=0.692, nonboundary=True
  omega: delta=0.6672, z=13.34; phi: delta=0.0490, z=0.29; beta: delta=0.2900, z=29.00
- QW-1741: final_score=0.0000, raw_score=0.0000, fit=0.752, nonboundary=False
  omega: delta=0.6354, z=21.18; phi: delta=1.1236, z=7.33; beta: delta=0.2400, z=24.00
- QW-1743: final_score=0.0000, raw_score=0.0000, fit=0.868, nonboundary=False
  omega: delta=0.5854, z=19.51; phi: delta=0.7461, z=9.33; beta: delta=0.1773, z=17.73

## Identifiability Penalty
- QW-1742 condition number: 3.431e+12
- QW-1744 condition number: 2.970e+12
- ident_penalty: 0.352

## Summary
- mean_base_score: 0.0000
- compatibility_index: 0.0000
- strongest_source: QW-1739 (0.0000)

## Artifacts
- JSON: `report_qw1862_micromodel_canonical_compatibility_audit.json`
