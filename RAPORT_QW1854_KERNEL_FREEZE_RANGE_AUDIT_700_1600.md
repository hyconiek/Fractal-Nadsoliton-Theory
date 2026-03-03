# RAPORT QW-1854: KERNEL FREEZE RANGE AUDIT (QW-700..1600)

- Data UTC: 2026-03-02T23:44:16.563002+00:00
- Verdict: **RANGE_700_1600_NOT_FULLY_TRACED**

## Scan Stats
- files_scanned: 7556
- files_with_qw_range_refs: 442
- files_kernel_related_in_range: 232

## Coverage
- any-ref coverage: 537/901 (0.596)
- kernel-ref coverage: 534/901 (0.593)

## Contradictions
- phi dual definition (pi/6 vs 0): True
- beta dual definition (0.01 vs 0.05): True
- node-set conflict (2,5,8,11 vs 2,8,14): True
- contradiction_count: 3

## Notes
- This is an evidence-range audit, not a proof of each individual QW computation.

## Artifacts
- JSON: `report_qw1854_kernel_freeze_range_audit_700_1600.json`
