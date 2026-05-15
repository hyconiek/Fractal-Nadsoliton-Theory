# P1784 — S734
## STRICT CURRENT PRIORITY BLOCKER -> ACTION CONTRACT (PL)

## Cel

Przekształcić status blokad priorytetowych w jednoznaczny kontrakt działań
(`trigger -> delivery -> success_flag`) bez skrótów i bez fałszywych PASS.

## Technical progress

- Zmapowano bieżące blokery na trzy akcje wykonawcze (`A1`, `A2`, `A3`).
- Każda akcja ma jawny trigger i warunek sukcesu.
- Zachowano single-lane porządek zgodny z wcześniejszym freeze/exit matrix.

## Co zostało dowiedzione

1. Blokady nie są już tylko opisem statusu — mają przypisane kontrakty wykonawcze.
2. Kolejność `A1 -> A2 -> A3` jest jawna i zgodna ze strict no-false-pass.

## Co nadal jest OPEN

1. Domknięcie W1 FULL_EXPORT.
2. Joint componentwise witness H1+metric.
3. Formalne odblokowanie BW/BRST/CUT.

## Ryzyka false-pass

- Próba wykonania `A3` bez zamknięcia `A1` i `A2`.
- Deklaracja sukcesu akcji bez jawnego witness trace.

## Następny uczciwy krok

Wykonać `A1`, następnie `A2`, następnie `A3`, aktualizując status tylko na podstawie
udokumentowanych witnessów i polityki `PASS_ZERO/OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Wyjaśnienie dla laika

To zamiana „listy problemów” na „plan naprawczy”: dla każdego problemu jest
konkretny następny ruch i jasny sygnał, kiedy można uznać go za rozwiązany.
