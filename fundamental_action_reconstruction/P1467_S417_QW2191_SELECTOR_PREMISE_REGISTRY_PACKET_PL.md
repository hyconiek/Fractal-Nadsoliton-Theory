# P1467 — S4.17 QW-2191 Selector Premise Registry (PL)

Status: `P1467_EXECUTED_QW2191_PREMISE_REGISTRY_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

W odpowiedzi na zainteresowanie rozwiązaniem `QW-2191` przygotować uczciwy,
jawny rejestr kandydatów premis selektorowych bez udawania strict-core closure.

## Zasada rygoru

- bez legacy bridge,
- bez twierdzenia strict-core closure,
- każdy kandydat premisy ma status `NON_STRICT_UNLESS_PROVEN_INTERNAL`.

## Wynik docelowy

Eksport tabeli kandydatów i decyzji, który kandydat jest najlepszy do następnego
lokalnego testu (tylko jako proposal, nie jako domknięcie QW-2191).
