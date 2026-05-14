# P1604 / S554 — Noncyclic provider upgrade for G1 + bridge

## Cel
Po nieudanym skanie adaptacyjnym (P1603) wprowadzić plan nowej klasy providera noncyclic,
zgodnej z rygorem strict-only i QW-238x, aby odblokować G1 oraz bridge theorem.

## Wejścia
- `generated/p1603_s553_adaptive_tail_depth_and_bridge_cotuning_summary.json`

## Wyjście
- `generated/p1604_s554_noncyclic_provider_upgrade_for_g1_bridge_summary.json`

## Rygor
- Bez legacy bridge.
- Bez walidacji zewnętrznej.
- Jawne wskazanie nowego blocker-cut oraz noncyclic anchor.
