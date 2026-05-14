# P1601 / S551 — Strict selector-source tail correction candidate

## Cel
Wykonać wewnętrzny (strict-only) kandydat korekty ogona dla najgorszych sektorów G1,
aby sprawdzić czy można odblokować witness `W_G1_full_domain_selector_gap_discharge`.

## Wejścia
- `generated/p1600_s550_targeted_g1_blocker_sectorization_summary.json`

## Wyjście
- `generated/p1601_s551_strict_selector_source_tail_correction_candidate_summary.json`

## Rygor
- Bez legacy bridge.
- Bez walidacji zewnętrznej.
- Jawne porównanie metryk przed/po korekcie kandydującej.
