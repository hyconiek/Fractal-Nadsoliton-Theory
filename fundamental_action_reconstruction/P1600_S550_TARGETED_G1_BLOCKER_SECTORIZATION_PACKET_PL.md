# P1600 / S550 — Targeted G1 blocker sectorization

## Cel
Zidentyfikować najgorsze sektory domeny, które blokują witness `W_G1_full_domain_selector_gap_discharge`,
aby następny krok był celowaną korektą strict selector-source (a nie powtórką całego pipeline).

## Wejścia
- `generated/p1581_s531_strict_selector_source_samples.csv`
- `generated/p1599_s549_replay_g3_with_refined_g1_witness_summary.json`

## Wyjście
- `generated/p1600_s550_targeted_g1_blocker_sectorization_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Bez walidacji zewnętrznej.
- Raportujemy jawnie najgorsze pary sektorów i metrykę tail readiness.
