# P1588 / S538 G1 Full-Domain Selector Gap Discharge Packet (PL)

Status: `P1588_EXECUTED_G1_FULL_DOMAIN_CHECK`
As of: `2026-05-14`

## Cel

Po `P1587` zaatakować `G1` na pełnej domenie strict:

1. policzyć błąd selector-gap na całej domenie,
2. sklasyfikować pass/fail dla `G1`,
3. utrzymać tor strict-only bez bridge do legacy.

## Wynik

- Eksportuje metryki `full_domain_error_max` i `full_domain_error_l2`.
- Wskazuje, czy można przejść do `G2`, czy trzeba dalej rafinować `G1`.

## Artefakt

- `generated/p1588_s538_g1_full_domain_selector_gap_discharge_summary.json`

## Następny uczciwy krok

`P1589`: złożyć `G1` z `G2` (global stability) albo zawęzić klasy kontrprzykładów dla `G1`.
