# P1528 S478 Theorem-Level Selector Uniqueness Scaffold Packet (No Legacy Bridge)

Status: `P1528_EXECUTED_THEOREM_LEVEL_SELECTOR_UNIQUENESS_SCAFFOLD_PROVISIONAL`
As of: `2026-05-13`

## Cel

Następny uczciwy krok po `P1527`:

- przejść z heurystyki wyboru gałęzi do formalnego szkicu dowodowego,
- utrzymać strict-only,
- nie zamykać `QW-2191` bez pełnego dowodu.

## Zakres

`S478` buduje szkic twierdzenia (scaffold), który weryfikuje trzy warunki:

1. **Existence**: istnieje gałąź kandydująca po redukcji,
2. **Dominance**: gałąź wybrana dominuje punktowo nad alternatywami w jawnej
   metryce,
3. **Stability**: dominacja utrzymuje się dla zestawu perturbacji wejściowych.

## Kontrakt wyjścia

- `existence_pass` (bool),
- `dominance_pass` (bool),
- `stability_pass` (bool),
- `theorem_scaffold_pass` (bool = koniunkcja trzech warunków),
- `qw2191_closed` (pozostaje `false` na etapie scaffold).

## PASS/FAIL

PASS jeśli:

1. wszystkie trzy warunki są jawnie liczone,
2. metryka dominacji jest reprodukowalna,
3. scaffold nie jest mylony z pełnym twierdzeniem.

FAIL jeśli:

1. którykolwiek warunek jest implicit,
2. wynik scaffold automatycznie zamyka `QW-2191`,
3. pojawia się legacy transfer.
