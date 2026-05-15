# P1721 S671 Strict Componentwise Curvature Variation Pack Packet (PL)

Status: `P1721_EXECUTED_STRICT_CURVATURE_VARIATION_PACK_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1720` dostarczyć brakujący pakiet wariacji krzywizny potrzebny do
wykonania metrycznego testu residual w runnerze componentwise.

## Co wyeksportowano

1. Wzory `δ(R^2)`, `δ(Ricci^2)`, `δ(Riemann^2)`.
2. Wzory pomocnicze `δR`, `δR_{ab}`, `δΓ`.
3. Notę integracyjną do runnera `componentwise_sympy_metric_residual_runner_v1`.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To dokładnie te brakujące „klocki matematyczne”, które wcześniej blokowały
obliczenie równania grawitacyjnego. Teraz można zrobić kolejny krok i policzyć
rzeczywistą resztę testu.

## Następny uczciwy krok (rekomendacja)

Wpiąć ten pakiet do runnera i opublikować pierwszy wynik `ELg-Emunu`:
`zero` albo `obstrukcja`.
