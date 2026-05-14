# P1530 S480 Strict-Core External Reproduction Packet (No Legacy Bridge)

Status: `P1530_EXECUTED_STRICT_CORE_EXTERNAL_REPRODUCTION_PACKET_PROVISIONAL`
As of: `2026-05-13`

## Cel

Następny uczciwy krok po `P1529`:

- przejść z reprodukcji wewnętrznej do reprodukcji zewnętrznej,
- sprawdzić zgodność wyniku scaffold dla wejścia dostarczonego spoza bieżącego
  runu,
- utrzymać strict-only, bez legacy bridge.

## Zakres

`S480` eksportuje minimalny protokół zewnętrznej reprodukcji:

1. przyjęcie `external_case_id` i `external_branch_set`,
2. uruchomienie tych samych testów dominance/stability,
3. porównanie z runem referencyjnym,
4. jawne oznaczenie `qw2191_closed=false`.

## Kontrakt wyjścia

- `reference_run`,
- `external_run`,
- `cross_run_margin_delta`,
- `external_reproduction_pass`,
- `qw2191_closed=false`.

## PASS/FAIL

PASS jeśli zewnętrzny run przechodzi scaffold i jest zgodny z referencją
w granicach tolerancji.

FAIL jeśli reprodukcja zewnętrzna nie przechodzi albo jeśli wynik jest mylony
z pełnym domknięciem selektora.
