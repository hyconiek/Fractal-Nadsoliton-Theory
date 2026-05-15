# P1773 / S723 — Synchronizacja state-vector po planie W1–W4

Status: `P1773_S723_REVERSE_GATE_STATE_SYNC_WITH_PROMOTION_GUARD_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Po `P1772` wykonano synchronizację state-vectora reverse-gates, aby formalnie
związać stan bramek z planem W1–W4.

Nowy stan:

- `G_BW`: aktywna obstrukcja, wymaga W1–W4,
- `G_BRST`: zablokowane przez `G_BW`,
- `G_CUT`: zablokowane przez `G_BW` i `G_BRST`,
- `qg_theorem_promotable`: `False`.

## Co zostało dowiedzione

1. Mamy spójny i jawny guard promocji: brak możliwości theorem-level update,
   dopóki `G_BW != PASS_ZERO` i W1–W4 nie są dowiezione jako `FULL_EXPORT`.
2. State-vector jest teraz zsynchronizowany z realnym bottleneckiem,
   a nie z ogólną deklaracją gotowości.

## Co nadal jest OPEN

1. W1–W4 (pełne eksporty tensorowe) nadal otwarte.
2. Kolejne wykonanie `G_BW`.
3. Odblokowanie BRST/Cutkosky i dalsze bramki QG.

## Ryzyka false-pass

1. "Papierowe" podniesienie readiness bez domknięcia W1–W4.
2. Pominięcie rozróżnień `LOCAL/GLOBAL`, `REDUCED/NONPROXY`, `SCAFFOLD/FULL_EXPORT`.
3. Promocja BRST/Cutkosky przed finalnym werdyktem `G_BW`.

## Następny uczciwy krok

Dowiezienie W1–W4 i uruchomienie kolejnej próby `G_BW` z pełnym trace.
Dopiero po realnym `PASS_ZERO` aktualizacja BRST/Cutkosky.

## Dla laika

To aktualizacja "tablicy kontrolnej": system jasno pokazuje, że są konkretne
prace do ukończenia i nie da się przeskoczyć do finałowych testów skrótem.
