# P1771 / S721 — Retry G_BW po sanitizacji notacji

Status: `P1771_S721_GBW_RETRY_AFTER_SANITIZATION_OBSTRUCTION_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Wykonano retry `G_BW` po sanitizacji notacji z `P1770`.

Wynik retry:

`OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

Różnica względem `P1769`:

- błąd notacyjno-parserowy został usunięty,
- obstrukcja pozostała, ale jest teraz jednoznacznie matematyczna (nie parserowa).

## Co zostało dowiedzione

1. Sanitizacja notacji była konieczna i skuteczna technicznie.
2. Po usunięciu problemu notacji główna blokada nadal leży w rozwinięciach
   komponentowych krzywizny i kontrterminów.
3. Brak podstaw do `PASS_ZERO`; werdykt uczciwie pozostaje `OBSTRUCTION`.

## Co nadal jest OPEN

1. Pełne rozwinięcia `H_R2/H_Ric2/H_Riem2`.
2. Pełny komponentowy basis `T_CT^{mu nu}`.
3. Zamknięcie residualu `B1/B2/B3/C1/C2` i divergence trace.
4. Odblokowanie `G_BRST` i `G_CUT`.

## Ryzyka false-pass

1. Traktowanie „retry done” jako postępu theorem-level.
2. Lokalna poprawa techniczna mylona z fizycznym domknięciem.
3. Próba otwarcia BRST/Cutkosky mimo aktywnej obstrukcji `G_BW`.

## Następny uczciwy krok

Dowieźć pełny komponentowy basis tensorowy (`H_R2/H_Ric2/H_Riem2` + `T_CT`) i
powtórzyć `G_BW` na tej samej rodzinie teł oraz tej samej bazie
`B1/B2/B3/C1/C2`, bez zmiany kryterium werdyktu.

## Dla laika

To jak powtórzenie testu po naprawieniu błędu w oprogramowaniu: teraz wiemy,
że problem nie wynika z zapisu danych, tylko z trudnej części obliczeń,
którą trzeba jeszcze dopracować.
