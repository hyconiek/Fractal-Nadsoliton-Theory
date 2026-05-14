# P1603 / S553 — Adaptive tail depth and bridge co-tuning scan

## Cel
Przeskanować adaptacyjnie głębokość i siłę strict-tail correction dla G1,
oraz równocześnie oznaczyć konieczność współstrojenia theorem-bridge przed finalnym G3.

## Wejścia
- `generated/p1581_s531_strict_selector_source_samples.csv`
- `generated/p1602_s552_apply_tail_correction_and_replay_summary.json`

## Wyjście
- `generated/p1603_s553_adaptive_tail_depth_and_bridge_cotuning_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Bez walidacji zewnętrznej.
- Jawny skan siatki parametrów i raport kandydatów admissible.
