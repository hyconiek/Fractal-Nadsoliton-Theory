# RAPORT QW-1831: EVENT-WINDOWED GW FEATURES

- Data UTC: 2026-03-02T22:14:08.960208+00:00
- GPS: 1266965117
- n samples: 2097152
- Window/step [s]: 8.0/2.0
- Band [Hz]: 20.0-400.0
- Verdict: **EVENT_WINDOW_FEATURE_BASELINE_PARTIAL**

## Pair Summary (selected)
- H1-L1: n=253, median corr@10ms=-0.000371, median max|corr|=0.008429, P(corr10>0.02)=0.012
- H1-V1: n=253, median corr@10ms=-0.000027, median max|corr|=0.001444, P(corr10>0.02)=0.000
- L1-V1: n=253, median corr@10ms=0.000823, median max|corr|=0.006254, P(corr10>0.02)=0.000

## Consistency
- std pair median corr@10ms: 0.000502
- std pair median max|corr|: 0.002919

## Pass Flags
- enough_windows: True
- lag_signal_presence_h1l1: False
- inter_pair_consistency: True

## Artifacts
- JSON: `report_qw1831_event_windowed_gw_features.json`
- CSV: `gw1831_window_features.csv`
