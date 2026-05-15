# P1719 S669 Strict Componentwise SymPy Runner Scaffold Packet (PL)

Status: `P1719_EXECUTED_STRICT_COMPONENTWISE_RUNNER_SCAFFOLD_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1718` przygotować jawny scaffold runnera `componentwise_sympy` do testu
metrycznego `EL_g - E_munu`.

## Co wyeksportowano

1. Kontrakt wejść/wyjść runnera.
2. Kolejność kroków obliczeniowych.
3. Status gotowości do pierwszego uruchomienia.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To „instrukcja działania” dla narzędzia, które ma policzyć najtrudniejsze równanie
w teorii. Dzięki temu następny krok może już dać konkretny wynik zamiast planu.

## Następny uczciwy krok (rekomendacja)

Uruchomić runner na pierwszym ansatzie metryki i opublikować wynik residual:
`zero` albo `obstrukcja`.
