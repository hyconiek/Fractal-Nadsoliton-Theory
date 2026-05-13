# P1376 Full c_i Export For Coefficient Class (PL)

Status: `P1376_EXECUTED_CLASS_SCHEMA_EXPORT_PARTIAL_NO_FALSE_PASS`
As of: `2026-05-12`

## Cel

Po `P1375` wykonujemy krok konieczny do realnego testu `L-B1-02`:

jawny eksport klasy współczynników `c_i` używanych w
`transport_drift = max_i |c_i(A)-c_i(B)| / max(1,|c_i(A)|)`.

## Zakres

Ten pakiet nie domyka jeszcze inwariancji transportu.
Dostarcza wyłącznie formalną klasę i schemat danych `c_i`.

## Minimalna klasa v1

Eksportujemy klasę indeksów:

- `c_g3`  (slot SU(3)),
- `c_g2`  (slot SU(2)),
- `c_g1`  (slot U(1)),
- `c_mix` (slot mieszania gauge).

## Status semantyczny

- schema class: `DEFINED`,
- complete numerical population on strict lane: `NOT_YET`,
- theorem-level invariance result: `OPEN`.

## Decyzja profesorska

To jest wymagany krok infrastrukturalny, ale nie wolno mieszać go
z dowodem fizycznym. Dowód `L-B1-02` zaczyna się dopiero po
pełnym wypełnieniu tej klasy rzeczywistymi danymi z pipeline.

