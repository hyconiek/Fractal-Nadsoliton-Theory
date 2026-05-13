# P1492 — S4.42 QW-2191 Selector Theorem Candidate (PL)

Status: `P1492_EXECUTED_QW2191_SELECTOR_THEOREM_CANDIDATE_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Po robust sweep (`P1491`) przygotować pierwszy kandydat twierdzenia
strict-core o źródle selektora, bez legacy bridge.

## Decyzja profesorska

Budujemy kandydat twierdzenia tylko na robust podzakresie `kappa` i przy
jawnych założeniach:

1. `|Delta_SB(kappa)| <= safety_margin`,
2. `G1(kappa) < G0` na całym podzakresie,
3. brak zmiany znaku orientacji selektora,
4. brak ukrytych importów legacy.

To nadal nie jest final closure theorem, ale formalny szkic dowodowy
z warunkami falsyfikowalnymi.
