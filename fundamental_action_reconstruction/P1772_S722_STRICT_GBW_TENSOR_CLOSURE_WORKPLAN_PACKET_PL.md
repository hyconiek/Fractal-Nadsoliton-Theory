# P1772 / S722 — Plan domknięcia tensorowego dla G_BW

Status: `P1772_S722_GBW_TENSOR_CLOSURE_WORKPLAN_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Po `P1771` (retry nadal z obstrukcją) ustalono precyzyjny plan domknięcia
czterech brakujących bloków tensorowych, bez zmiany tła i bez zmiany polityki
werdyktu.

Pakiet rozbija blocker na zadania W1–W4:

1. `W1`: pełny komponentowy export `H_R2`,
2. `W2`: pełny komponentowy export `H_Ric2` z komutatorami pochodnych,
3. `W3`: pełny komponentowy export `H_Riem2` z transportem podwójnej dywergencji,
4. `W4`: pełny komponentowy basis `T_CT` zgodny z `B1/B2/B3/C1/C2`.

## Co zostało dowiedzione

1. Obstrukcja `G_BW` ma już dokładną dekompozycję roboczą (W1-W4),
   więc praca może iść równomiernie i audytowalnie.
2. Kontrakt re-execution jest zamrożony:
   - stała rodzina teł,
   - stała baza residualu,
   - stała polityka werdyktu (`PASS_ZERO`/`OBSTRUCTION_WITH_DIVERGENCE_TRACE`).

## Co nadal jest OPEN

1. Wykonanie i dostarczenie W1-W4.
2. Kolejne wykonanie `G_BW` po W1-W4.
3. Odblokowanie `G_BRST` i `G_CUT`.

## Ryzyka false-pass

1. Częściowe domknięcie jednego bloku (np. tylko W1) nie może być mylone z
   pełnym domknięciem `G_BW`.
2. Nie wolno zmieniać klasy tła/bazy przy kolejnym uruchomieniu.
3. Nie wolno promować BRST/Cutkosky przed finalnym werdyktem `G_BW`.

## Następny uczciwy krok

Dostarczyć W1-W4 i uruchomić kolejne `G_BW` z pełnym residual/divergence trace,
bez zmiany kontraktu decyzji.

## Dla laika

Zamiast zgadywać, rozpisaliśmy problem na cztery konkretne zadania matematyczne.
Gdy każde z nich zostanie ukończone, test spójności będzie można powtórzyć w
sposób uczciwy i porównywalny.
