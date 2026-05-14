# P1598 / S548 — Strict G1 witness refinement (full-domain selector gap)

## Cel
Wzmocnić brakujący witness `W_G1_full_domain_selector_gap_discharge` surowszym testem globalnym,
pozostając całkowicie w torze strict-only.

## Wejścia
- `generated/p1581_s531_strict_selector_source_samples.csv`
- `generated/p1597_s547_final_g3_theorem_composition_object_summary.json`

## Wyjście
- `generated/p1598_s548_strict_g1_witness_refinement_summary.json`

## Rygor
- Bez bridge do legacy.
- Bez walidacji zewnętrznej.
- PASS tylko przy spełnieniu jednoczesnych progów dla `max`, `l2` i `mean` błędu full-domain.
