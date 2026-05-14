# P1658 / S608 — Strict H1 package integration and H2 bootstrap

## Cel
Scalić witnessy H1 (`phi-H`, `gauge-scalar`, `gauge-metric`) w jeden pakiet operatorowy
oraz wykonać bootstrap H2 dla tego samego bloku, bez false-pass.

## Zakres
- strict-only, bez bridge do legacy,
- pełny tor kontekstowy: `K_strict -> współczynniki -> L_total -> EOM`,
- tor odwrotny pozostaje `OPEN` aż do domknięcia H1..H4 i bram QG.

## Konstrukcja
- Integracja H1:
  - P1655: `phi-H` local,
  - P1656: `gauge-scalar` covariant,
  - P1657: `gauge-metric` covariant.
- Bootstrap H2:
  - zapis warunku self-adjointness dla liniaryzacji operatora EOM,
  - status: `PARTIAL` (schema + obligations).

## Wyjście
- `generated/p1658_s608_strict_h1_package_and_h2_bootstrap_summary.json`
