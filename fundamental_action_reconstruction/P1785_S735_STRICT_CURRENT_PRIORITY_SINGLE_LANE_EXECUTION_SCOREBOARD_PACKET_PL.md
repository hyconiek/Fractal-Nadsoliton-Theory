# P1785 — S735
## STRICT CURRENT PRIORITY SINGLE-LANE EXECUTION SCOREBOARD (PL)

## Cel

Opublikować prostą tablicę wykonawczą: który krok jest wykonywalny teraz,
a które są zablokowane przez zależności upstream.

## Technical progress

- Wyeksportowano scoreboard `L1 -> L2 -> L3`.
- Ustalono flagi `ready_to_execute` dla każdego kroku.
- Zachowano single-lane porządek zgodny z kontraktem blokad.

## Co zostało dowiedzione

1. Tylko `L1/A1` jest wykonywalny natychmiast.
2. `L2` i `L3` pozostają zablokowane do czasu spełnienia warunków upstream.

## Co nadal jest OPEN

1. W1 FULL_EXPORT (L1 completion condition).
2. Joint H1+metric witness run (L2 completion condition).
3. Release theorem-gate freeze (L3 completion condition).

## Ryzyka false-pass

- Ominięcie `L1` i przejście bezpośrednio do `L2/L3`.
- Oznaczenie kroku jako zakończony bez spełnienia jawnych completion conditions.

## Następny uczciwy krok

Wykonać `L1`, udokumentować warunek domknięcia, następnie formalnie odblokować `L2`.

## Wyjaśnienie dla laika

To harmonogram z blokadami: najpierw pierwszy krok, dopiero potem drugi i trzeci.
