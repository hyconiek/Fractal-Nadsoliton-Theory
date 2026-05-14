# P1531 S481 Strict-Core Peer-Review Reproduction Packet (No Legacy Bridge)

Status: `P1531_EXECUTED_STRICT_CORE_PEER_REVIEW_REPRODUCTION_PROVISIONAL`
As of: `2026-05-14`

## Cel

Następny uczciwy krok po `P1530`:

- przejść od pojedynczej reprodukcji zewnętrznej do wielo-case reprodukcji
  peer-review,
- policzyć rozrzut metryk między niezależnymi case'ami,
- utrzymać strict-only i brak legacy bridge.

## Zakres

`S481` uruchamia min. 3 zewnętrzne case'y i raportuje:

1. `cross_run_margin_delta` dla każdego case'u,
2. `delta_mean` i `delta_max`,
3. `peer_review_reproduction_pass` przy jawnej tolerancji.

## Kontrakt wyjścia

- `reference_run`,
- `external_cases_results`,
- `delta_mean`,
- `delta_max`,
- `peer_review_reproduction_pass`,
- `qw2191_closed=false`.

## PASS/FAIL

PASS jeśli wszystkie case'y przechodzą scaffold i `delta_max <= tolerance`.

FAIL jeśli dowolny case nie przechodzi lub jeśli wynik jest mylony z pełnym
strict selector closure.
