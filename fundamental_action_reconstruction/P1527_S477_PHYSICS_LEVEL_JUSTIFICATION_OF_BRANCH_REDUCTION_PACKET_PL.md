# P1527 S477 Physics-Level Justification of Branch Reduction Packet (No Legacy Bridge)

Status: `P1527_EXECUTED_PHYSICS_LEVEL_JUSTIFICATION_OF_S476_RULE_PROVISIONAL`
As of: `2026-05-13`

## Cel

Następny uczciwy krok po `P1526`:

- podnieść regułę redukcji gałęzi z poziomu "czysto techniczny filtr"
  do poziomu "fizycznie uzasadniona heurystyka strict-only",
- bez bridge do legacy,
- bez domknięcia `QW-2191`.

## Zakres S477

`S477` nie dowodzi pełnej teorii selektora. Dostarcza tylko:

1. jawny scoring fizyczny gałęzi (`branch_physics_score`),
2. deterministyczny wybór gałęzi maksymalnie zgodnej z rygorem strict-only,
3. test stabilności wyboru na wielu zestawach wejściowych.

## Minimalny scoring fizyczny (heurystyczny)

Dla każdej gałęzi obliczamy:

```text
score = w_anchor*I(strict_anchor_) + w_noncyclic*I(noncyclic_marker) + w_trace*I(trace_marker)
```

z wagami stałymi i jawnymi. To pozostaje heurystyka, ale już audytowalna.

## PASS/FAIL

PASS jeśli:

1. wybór gałęzi jest reprodukowalny,
2. scoring i wagi są jawne,
3. stabilność wyboru utrzymuje się na zadanych zestawach,
4. `qw2191_closed` pozostaje `false` (brak pełnego theorem-level witness).

FAIL jeśli:

1. wybór jest niereprodukowalny,
2. reguła zależy od ukrytego parametru,
3. pojawia się claim strict closure bez dodatkowego twierdzenia.
