# RAPORT QW-1799: HIERARCHICAL TRANSFER TEST

- Data UTC: 2026-03-02T20:49:28.185019+00:00
- Family: mixed_low [1, 2, 3]
- Shrink factor: 0.20
- Splits: 14
- Test P(hier>M2): 0.857
- Test P(hier>flat): 0.786
- Test median delta hier-M2: 1.0168
- Test std delta hier-M2: 1.1192
- Gap median (train-test): 0.1083
- Verdict: **HIERARCHICAL_TRANSFER_WEAK**

## Pass Flags
- transfer_gain: True
- transfer_stability: False
- generalization_gap_control: True
- test_positive_vs_flat: False

## Artifacts
- JSON: `report_qw1799_hierarchical_transfer_test.json`
