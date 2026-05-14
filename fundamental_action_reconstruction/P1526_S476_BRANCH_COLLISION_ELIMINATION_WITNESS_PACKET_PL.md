# P1526 S476 Branch-Collision Elimination Witness Packet (No Legacy Bridge)

Status: `P1526_EXECUTED_BRANCH_COLLISION_ELIMINATION_WITNESS_NO_FALSE_CLOSURE`
As of: `2026-05-13`

## Cel

Następny uczciwy krok po `P1525`: dodać minimalny świadek eliminacji kolizji
gałęzi selektora strict-core, aby test unikalności nie operował na surowej
wielo-gałęziowości.

Trasa:

```text
F_nadsoliton -> L_SM + L_GR
```

## Założenie rygoru

1. Brak legacy bridge.
2. Brak automatycznego domknięcia `QW-2191`.
3. Najpierw jawna redukcja gałęzi (`branch reduction`), potem test unikalności.

## Minimalny kontrakt

Wejście:

- `alternative_branch_ids`,
- `reduction_rule_id`.

Wyjście:

- `reduced_branch_ids`,
- `elimination_applied` (bool),
- `elimination_reason_code`.

## Reguła minimalna S476

Dla tego checkpointu stosujemy konserwatywną regułę:

- jeśli istnieje gałąź oznaczona prefiksem `strict_anchor_`, wybieramy ją jako
  jedyną gałąź zredukowaną,
- w przeciwnym razie nie eliminujemy gałęzi.

To jest krok techniczny, nie dowód pełnej fizycznej wystarczalności.

## PASS/FAIL

PASS jeśli:

1. redukcja jest audytowalna (`reduction_rule_id` + reason code),
2. wynik unikalności liczony jest na `reduced_branch_ids`,
3. `qw2191_closed` pozostaje zależne od `intake_pass AND uniqueness_pass`.

FAIL jeśli:

1. redukcja jest niejawna,
2. unikalność liczona jest na nieprzefiltrowanych gałęziach mimo aktywnej
   reguły redukcji,
3. pojawia się closure claim bez przejścia wszystkich bramek.
