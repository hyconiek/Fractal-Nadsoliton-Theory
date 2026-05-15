# P1793 S743 Strict SV1->SV5 Evidence Intake Validator Packet (PL)

## Cel

Zrealizować kolejny krok po `P1792`: dostarczyć działający walidator intake,
który automatycznie ocenia poprawność schemy i spójność freeze dla ledgerów `SV1..SV5`.

## Zakres

Walidator przyjmuje JSON z listą ledgerów i wykonuje:

1. walidację pól wymaganych,
2. walidację `result_kind`,
3. walidację jawnego `residual_vector` przy `PASS_ZERO`,
4. walidację spójności `background/index/boundary` między wszystkimi ledgerami,
5. zwrot deterministycznego werdyktu intake.

## Werdykt intake

- `PASS_ZERO` tylko gdy:
  - wszystkie ledgery są poprawne schema,
  - freeze jest spójny,
  - brak naruszeń `INVALID_PASS_CLAIM`.

- W przeciwnym razie:
  - `OPEN_OBSTRUCTION_WITH_TRACE` i jawna lista tagów obstrukcji.

## Obstrukcje

- `INVALID_LEDGER_SCHEMA`
- `INVALID_RESULT_KIND`
- `INVALID_PASS_CLAIM`
- `FREEZE_MISMATCH`

## Zasada no-false-pass

Ten walidator nie zamyka theorem gates. To wyłącznie kontrola jakości wejścia
przed protokołem `P1791`.

## Następny uczciwy krok

Podłączyć realne ledgery wykonania `SV1..SV5` do walidatora i publikować
wynik intake razem z obstruction trace.

## Objaśnienie dla laika

To automat sprawdzający „czy raporty są poprawnie wypełnione i dotyczą tych samych warunków testu”, zanim w ogóle zaczniemy oceniać wynik fizyczny.
