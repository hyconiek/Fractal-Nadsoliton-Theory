# P1529 S479 Independent Strict-Core Reproduction Packet (No Legacy Bridge)

Status: `P1529_EXECUTED_INDEPENDENT_STRICT_CORE_REPRODUCTION_PROVISIONAL`
As of: `2026-05-13`

## Cel

Następny uczciwy krok po `P1528`:

- niezależnie odtworzyć wynik scaffold (`existence/dominance/stability`) na
  drugim zestawie wejściowym,
- porównać zgodność metryk między przebiegiem A i B,
- utrzymać strict-only i brak legacy bridge.

## Zakres

`S479` daje "independent reproduction check", nie pełne zamknięcie.

Warunki:

1. oba przebiegi muszą przejść `theorem_scaffold_pass`,
2. różnica `dominance_margin` między przebiegami musi być <= tolerancji,
3. `qw2191_closed` pozostaje `false` do czasu niezależnego theorem-level proof.

## Kontrakt wyjścia

- `run_a`,
- `run_b`,
- `margin_delta`,
- `reproduction_pass`,
- `qw2191_closed=false`.

## PASS/FAIL

PASS jeśli odtworzenie jest zgodne w granicach tolerancji i jawnie oznaczone
jako etap przed-domknięciowy.

FAIL jeśli pojawia się claim closure wyłącznie z reprodukcji scaffold.
