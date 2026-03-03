# RAPORT QW-1855: KERNEL FREEZE DEEPREAD QUEUE (QW-700..1600)

- Data UTC: 2026-03-02T23:44:18.425431+00:00
- Verdict: **DEEPREAD_QUEUE_PREPARED**

## Summary
- missing ids: 364
- ids with trace: 537
- priority counts: {'P0_MISSING': 364, 'P1_HIGH_RISK': 179, 'P2_MEDIUM_RISK': 33, 'P3_LOW_RISK': 325}
- canonical candidates count: 69

## Top Queue (first 25)
- QW-1148: P0_MISSING | risk=1.0
- QW-1149: P0_MISSING | risk=1.0
- QW-1150: P0_MISSING | risk=1.0
- QW-1151: P0_MISSING | risk=1.0
- QW-1154: P0_MISSING | risk=1.0
- QW-1155: P0_MISSING | risk=1.0
- QW-1156: P0_MISSING | risk=1.0
- QW-1157: P0_MISSING | risk=1.0
- QW-1161: P0_MISSING | risk=1.0
- QW-1162: P0_MISSING | risk=1.0
- QW-1163: P0_MISSING | risk=1.0
- QW-1164: P0_MISSING | risk=1.0
- QW-1165: P0_MISSING | risk=1.0
- QW-1166: P0_MISSING | risk=1.0
- QW-1167: P0_MISSING | risk=1.0
- QW-1168: P0_MISSING | risk=1.0
- QW-1169: P0_MISSING | risk=1.0
- QW-1170: P0_MISSING | risk=1.0
- QW-1171: P0_MISSING | risk=1.0
- QW-1172: P0_MISSING | risk=1.0
- QW-1173: P0_MISSING | risk=1.0
- QW-1174: P0_MISSING | risk=1.0
- QW-1175: P0_MISSING | risk=1.0
- QW-1176: P0_MISSING | risk=1.0
- QW-1177: P0_MISSING | risk=1.0

## Top Canonical Candidates (first 20)
- QW-1300: frozen_refs=6, kernel_refs=8, any_refs=16
- QW-986: frozen_refs=6, kernel_refs=7, any_refs=11
- QW-1202: frozen_refs=5, kernel_refs=5, any_refs=12
- QW-1204: frozen_refs=4, kernel_refs=6, any_refs=10
- QW-1206: frozen_refs=4, kernel_refs=6, any_refs=13
- QW-1213: frozen_refs=4, kernel_refs=6, any_refs=13
- QW-1210: frozen_refs=4, kernel_refs=3, any_refs=7
- QW-1097: frozen_refs=3, kernel_refs=6, any_refs=7
- QW-1106: frozen_refs=3, kernel_refs=6, any_refs=7
- QW-1116: frozen_refs=3, kernel_refs=6, any_refs=7
- QW-1159: frozen_refs=3, kernel_refs=6, any_refs=11
- QW-1400: frozen_refs=3, kernel_refs=6, any_refs=14
- QW-978: frozen_refs=3, kernel_refs=5, any_refs=6
- QW-992: frozen_refs=3, kernel_refs=5, any_refs=6
- QW-911: frozen_refs=3, kernel_refs=4, any_refs=5
- QW-980: frozen_refs=3, kernel_refs=4, any_refs=5
- QW-988: frozen_refs=3, kernel_refs=4, any_refs=5
- QW-1000: frozen_refs=3, kernel_refs=4, any_refs=6
- QW-1100: frozen_refs=3, kernel_refs=4, any_refs=5
- QW-1122: frozen_refs=3, kernel_refs=4, any_refs=7

## Artifacts
- JSON: `report_qw1855_kernel_freeze_deepread_queue_700_1600.json`
