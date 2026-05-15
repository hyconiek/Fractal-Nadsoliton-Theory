# P1725 S675 Strict Componentwise Metric Residual Basis Reduction Pass1 Packet (PL)

Status: `P1725_EXECUTED_STRICT_METRIC_BASIS_REDUCTION_PASS1_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1724` wykonać pierwszy pass redukcji bazy tensorowej dla grupy
pochodnych kowariantnych krzywizny.

## Co wyeksportowano

1. Zastosowane reguły przepisania.
2. Baza przed/po redukcji.
3. Status: `BASIS_NORMALIZED_PASS1`.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To krok „sprzątania matematyki”: porządkujemy trudne składniki do wspólnego
języka, żeby w kolejnym kroku policzyć finalną resztę równania grawitacyjnego.

## Następny uczciwy krok (rekomendacja)

Wykonać pass2 dla składników Riemanna i złożyć residual na jednej bazie.
