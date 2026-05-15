# P1737 S687 Strict Full-Chain Next Action Contract (H1 then Metric) Packet (PL)

Status: `P1737_EXECUTED_NEXT_ACTION_CONTRACT_STRICT_ONLY`  
As of: `2026-05-15`

## Cel

Zamknąć etap planowania i przejść do jednoznacznej kolejności wykonawczej:

1. najpierw `H1`,
2. potem `EL_g - E_{μν}`,
3. dopiero potem aktualizacja bramek QG.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- decyzje tylko `PASS_ZERO` lub `OBSTRUCTION`.

## Następny uczciwy krok (rekomendacja)

Natychmiast wykonać `Step_1` z `P1737` i opublikować wynik.
