# P1591 / S541 Internal Replay Or Gap Focus Packet (PL)

Status: `P1591_EXECUTED_STATUS_TO_WORK_PACKAGE_TRANSLATION`
As of: `2026-05-14`

## Cel

Przetłumaczyć status po `P1590` na operacyjny pakiet pracy strict-only:

1. `CLOSED` -> pakiet replikacji wewnętrznej,
2. `OPEN` -> pakiet domykania braków theorem-level,
3. bez legacy bridge, bez fałszywego domknięcia.

## Wynik

- Eksportuje deterministyczny `work_package` zależny od statusu.
- Utrzymuje pełny tor strict i jawny status closure.

## Artefakt

- `generated/p1591_s541_internal_replay_or_gap_focus_summary.json`

## Następny uczciwy krok

Jeśli `OPEN`: wykonać focused theorem discharge; jeśli `CLOSED`: wykonać replay i zamrozić pakiet.
