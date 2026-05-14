# P1534 S484 LEM4 Counterexample Scanner Packet (No Legacy Bridge)

Status: `P1534_EXECUTED_LEM4_COUNTEREXAMPLE_SCANNER_PROVISIONAL`
As of: `2026-05-14`

## Cel

Następny uczciwy krok po `P1533`:

- zaatakować najważniejszą lukę `LEM_4_NO_EQUIVALENT_ALTERNATIVE_BRANCH`,
- najpierw próbą falsyfikacji (counterexample-first discipline),
- bez bridge do legacy.

## Zakres

`S484` uruchamia skaner kontrprzykładów dla klasy alternatywnych gałęzi i
sprawdza, czy istnieje gałąź o równoważnym `dominance_margin` względem gałęzi
wybranej.

Jeśli skaner znajdzie równoważną gałąź => `LEM_4` nie może być podniesiony.

## Kontrakt wyjścia

- `selected_branch`,
- `alternative_branches`,
- `equivalent_alternative_found` (bool),
- `lem4_status_update` in `{open, partial}`,
- `qw2191_closed=false`.

## PASS/FAIL

PASS jeśli skaner działa jawnie i daje reprodukowalny wynik.

FAIL jeśli `LEM_4` jest podniesiony bez negatywnego wyniku skanera.
