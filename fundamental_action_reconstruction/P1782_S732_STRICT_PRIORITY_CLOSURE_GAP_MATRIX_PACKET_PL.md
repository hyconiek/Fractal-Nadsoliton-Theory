# P1782 — S732
## STRICT PRIORITY CLOSURE GAP MATRIX (PL)

## Cel

Zamknąć lukę raportową między „mamy eksport operatora” a „mamy domknięcie theorem-level”
w postaci jednej matrycy klas braków.

## Technical progress

- Wyeksportowano matrycę braków dla priorytetów: `E_Aμ`, `E_H`, `EL_g`, `H1`, `BW`, `BRST/Cutkosky`.
- Każdy wpis ma klasę semantyczną zgodną ze strict no-false-pass.
- Matryca jest spięta z `P1781` (readiness delta).

## Co zostało dowiedzione

1. Operator-level nonproxy eksporty są utrzymane jako „done @ operator level”.
2. Wszystkie theorem-level bramki pozostają formalnie otwarte/zablokowane.

## Co nadal jest OPEN

1. Componentwise H1 witness.
2. Componentwise `EL_g-E_{μν}` witness.
3. Theorem-level `BW -> BRST -> Cutkosky`.

## Ryzyka false-pass

- Przepisanie klasy „operator done” na „closure done”.
- Aktualizacja bramek bez jawnego residual witness.

## Następny uczciwy krok

Domknąć W1 i wykonać joint componentwise witness run; aktualizować statusy wyłącznie
na podstawie `PASS_ZERO` albo `OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Wyjaśnienie dla laika

To tabela kontrolna jakości: pokazuje nie tylko postęp, ale też dokładnie czego
jeszcze brakuje do uczciwego finału.
