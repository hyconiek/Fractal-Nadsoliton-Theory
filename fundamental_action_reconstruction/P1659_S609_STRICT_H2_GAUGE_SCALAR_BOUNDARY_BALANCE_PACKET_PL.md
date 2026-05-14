# P1659 / S609 — Strict H2 gauge-scalar boundary-balance export

## Cel
Wykonać pierwszy jawny rachunek H2 po P1658:
`<u, L_E v> - <L_E u, v>` dla bloku gauge-scalar,
z rozdzieleniem części brzegowej i części objętościowej.

## Zakres
- strict-only, bez legacy bridge,
- wynik `PARTIAL` (brak globalnego domknięcia),
- status ToE/QG pozostaje `OPEN`.

## Konstrukcja
- Liniaryzowany blok operatora `L_E` dla pola cechowania i skalaru.
- Bilans formalny:
  - część objętościowa (bulk skew part),
  - całkowita dywergencja (boundary current).
- Celem checkpointu jest jawny eksport schematu rachunkowego i listy braków theorem-level.

## Wyjście
- `generated/p1659_s609_strict_h2_gauge_scalar_boundary_balance_summary.json`
