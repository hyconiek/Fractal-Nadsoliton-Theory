# RAPORT QW-1856: KERNEL FREEZE CONFLICT LEDGER (QW-700..1600)

- Data UTC: 2026-03-02T23:44:20.969792+00:00
- Verdict: **CONFLICT_LEDGER_PREPARED**

## Summary
- QW ids with trace: 540
- conflict rows (frozen + contradiction): 179
- canonical rows (frozen + no contradiction): 69

## Top Conflict Rows (first 20)
- QW-1200: c=3, frozen=10, kernel=17
- QW-700: c=2, frozen=6, kernel=14
- QW-857: c=2, frozen=1, kernel=6
- QW-877: c=2, frozen=1, kernel=4
- QW-896: c=2, frozen=1, kernel=4
- QW-842: c=1, frozen=9, kernel=7
- QW-722: c=1, frozen=7, kernel=11
- QW-983: c=1, frozen=5, kernel=7
- QW-826: c=1, frozen=4, kernel=8
- QW-977: c=1, frozen=4, kernel=7
- QW-707: c=1, frozen=4, kernel=6
- QW-996: c=1, frozen=4, kernel=6
- QW-828: c=1, frozen=4, kernel=5
- QW-840: c=1, frozen=4, kernel=5
- QW-957: c=1, frozen=4, kernel=5
- QW-827: c=1, frozen=4, kernel=4
- QW-829: c=1, frozen=4, kernel=4
- QW-830: c=1, frozen=4, kernel=4
- QW-831: c=1, frozen=4, kernel=4
- QW-833: c=1, frozen=4, kernel=4

## Top Canonical Rows (first 20)
- QW-1300: frozen=6, kernel=8, refs_any=16
- QW-986: frozen=6, kernel=7, refs_any=11
- QW-1202: frozen=5, kernel=5, refs_any=12
- QW-1204: frozen=4, kernel=6, refs_any=10
- QW-1206: frozen=4, kernel=6, refs_any=13
- QW-1213: frozen=4, kernel=6, refs_any=13
- QW-1210: frozen=4, kernel=4, refs_any=8
- QW-1097: frozen=3, kernel=7, refs_any=8
- QW-1106: frozen=3, kernel=7, refs_any=8
- QW-1116: frozen=3, kernel=7, refs_any=8
- QW-1159: frozen=3, kernel=7, refs_any=12
- QW-1400: frozen=3, kernel=7, refs_any=15
- QW-978: frozen=3, kernel=6, refs_any=7
- QW-992: frozen=3, kernel=6, refs_any=7
- QW-911: frozen=3, kernel=5, refs_any=6
- QW-980: frozen=3, kernel=5, refs_any=6
- QW-988: frozen=3, kernel=5, refs_any=6
- QW-1000: frozen=3, kernel=5, refs_any=7
- QW-1100: frozen=3, kernel=5, refs_any=6
- QW-1122: frozen=3, kernel=5, refs_any=8

## Artifacts
- JSON: `report_qw1856_kernel_freeze_conflict_ledger_700_1600.json`
