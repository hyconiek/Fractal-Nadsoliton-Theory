# P1720 S670 Strict Componentwise Runner First Execution Attempt Packet (PL)

Status: `P1720_EXECUTED_STRICT_COMPONENTWISE_FIRST_ATTEMPT_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1719` wykonać pierwsze uruchomienie runnera componentwise dla sektora metrycznego.

## Co wyeksportowano

1. Raport pierwszego uruchomienia runnera.
2. Status: `PARTIAL_EXECUTION_OBSTRUCTION`.
3. Jawna lista blokujących wariacji krzywiznowych.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To pierwsza próba odpalenia „silnika” dla grawitacji. Silnik ruszył częściowo,
ale zatrzymał się na najbardziej złożonych wzorach krzywizny. Teraz wiemy dokładnie,
co dopisać, żeby uzyskać pełny wynik.

## Następny uczciwy krok (rekomendacja)

Dopisać jawne wzory wariacyjne dla `R^2`, `Ricci^2`, `Riemann^2` w wersji componentwise,
a następnie policzyć `ELg-Emunu`.
