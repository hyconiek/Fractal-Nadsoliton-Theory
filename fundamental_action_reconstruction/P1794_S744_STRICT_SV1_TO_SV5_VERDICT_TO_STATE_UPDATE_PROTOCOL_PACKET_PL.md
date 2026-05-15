# P1794 S744 Strict SV1->SV5 Verdict-To-State Update Protocol Packet (PL)

## Cel

Zamknąć lukę między walidacją intake (`P1793`) a aktualizacją state-vector:
jak dokładnie mapować werdykt `PASS_ZERO`/`OPEN_OBSTRUCTION_WITH_TRACE` na statusy `SV1..SV5`.

## Wejścia

- `P1789` (bazowy state-vector),
- `P1793` (wynik walidacji intake i obstruction tags),
- `P1788` (theorem locks).

## Reguła mapowania

1. Jeśli `P1793.verdict == PASS_ZERO`:
   - `SV1..SV5 -> PASS_STRICT_LOCAL_EVIDENCE_ACCEPTED`.
2. Jeśli `P1793.verdict == OPEN_OBSTRUCTION_WITH_TRACE`:
   - `SV1..SV5 -> OPEN_OBSTRUCTION_WITH_TRACE`.

## Reguła globalna (niezmienna)

Niezależnie od wyniku lokalnego:

- `SV6..SV8` pozostają bez automatycznej promocji,
- theorem-level gates nadal zależą od globalnych witnessów BW/BRST/Cutkosky.

## Wymóg jawności

State update musi publikować:

- źródłowy werdykt intake,
- obstruction tags,
- nowe statusy `SV1..SV5`,
- niezmienione `SV6..SV8`.

## Następny uczciwy krok

Uruchomić generator checkpointu, który automatycznie produkuje ten update
bez ręcznej edycji statusów.

## Objaśnienie dla laika

To automatyczna reguła „jak przetłumaczyć wynik kontroli raportu na tablicę statusów projektu” — bez dowolnych interpretacji.
