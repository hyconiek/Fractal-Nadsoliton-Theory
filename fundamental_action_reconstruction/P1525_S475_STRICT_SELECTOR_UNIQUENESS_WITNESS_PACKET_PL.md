# P1525 S475 Strict Selector Uniqueness Witness Packet (No Legacy Bridge)

Status: `P1525_EXECUTED_STRICT_SELECTOR_UNIQUENESS_WITNESS_GATE_NO_FALSE_CLOSURE`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok po `P1524`:

- nie tylko intake-pass,
- ale minimalny test świadka unikalności selektora strict-core.

Trasa pozostaje:

```text
F_nadsoliton -> L_SM + L_GR
```

## Rygor

1. Brak bridge do legacy.
2. Brak automatycznej implikacji `intake_pass => qw2191_closed`.
3. `QW-2191` zamykamy tylko przy pozytywnym wyniku testu unikalności.

## Minimalny kontrakt świadka unikalności

Wejście:

- kandydat, który przeszedł `G_selector_intake^(strict)`,
- lista alternatywnych gałęzi selekcji (`alternative_branch_ids`).

Wyjście:

- `uniqueness_status` in `{unique_selector_confirmed, non_unique_or_unproven}`,
- `reason_code`,
- `branch_collision_count`.

## PASS/FAIL

PASS jeśli:

1. pipeline rozróżnia `intake_pass` i `uniqueness_pass`,
2. `qw2191_closed=true` tylko gdy `uniqueness_status=unique_selector_confirmed`,
3. dla przypadku niejednoznacznego utrzymuje `qw2191_closed=false`.

FAIL jeśli:

1. samo intake-pass zamyka `QW-2191`,
2. brak jawnej metryki kolizji gałęzi,
3. pojawia się claim closure bez świadka unikalności.
