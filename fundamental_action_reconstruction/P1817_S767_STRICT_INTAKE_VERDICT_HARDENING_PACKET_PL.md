# P1817 S767 Strict Intake Verdict Hardening Packet (PL)

Status: `P1817_EXECUTED_STRICT_INTAKE_VERDICT_HARDENING_NO_FALSE_PASS`

## Technical progress

Uszczelniono logikę `P1793` tak, by walidator nigdy nie zwracał `PASS_ZERO`,
jeśli choć jeden ledger ma `result_kind != PASS_ZERO`.

## Co zostało dowiedzione

1. Poprzednia reguła była proceduralnie zbyt słaba: poprawny schemat mógł dać `PASS_ZERO` nawet dla `OPEN_OBSTRUCTION_WITH_TRACE`.
2. Po poprawce walidator rozróżnia:
   - poprawność struktury danych,
   - fizyczny werdykt residualowy.
3. Aktualny wynik intake dla demo-ledgera to jawnie `OPEN_OBSTRUCTION_WITH_TRACE`.

## Co nadal OPEN

- Realny globalny witness dla `TG1_BW` na unified nonproxy runpack.
- Weryfikacja BRST/CUT po rzeczywistym odblokowaniu `TG1`.

## Ryzyka false-pass

1. Interpretacja "schema-valid" jako "physics-pass".
2. Używanie demo intake jako podstawy do promotion gate.
3. Pominięcie wymogu globalnego trace dla residualu.

## Następny uczciwy krok

Podstawić realny ledger z execution run (`P1806`) i ponowić `P1793 -> P1816 -> P1794`,
akceptując tylko binarny wynik z jawnie udokumentowanym residual witness.

## Krótkie wyjaśnienie dla laika

To poprawka „bezpiecznika”: nawet jeśli formularz jest poprawnie wypełniony,
to nie znaczy jeszcze, że fizyka się zgadza. Teraz system już tego nie pomyli.
