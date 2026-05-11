# P1215 W4 Gate Mode Comparison Packet

Status: `P1215_EXECUTED_STRICT_VS_PROVISIONAL_GATE_COMPARISON_NO_FALSE_PASS`
As of: `2026-05-11`

## Goal

Następny uczciwy krok po `P1214`: jawnie porównać zachowanie ścieżki W4
w trybie `STRICT_CERTIFIED` i `PROVISIONAL_UNCERTIFIED` na tych samych wejściach.

## Professor-level decision

Dodaję `p1215_w4_gate_mode_comparison.py`, który:

1. uruchamia `P1209` w obu trybach,
2. uruchamia `P1210` na wynikach obu bramek,
3. zapisuje jedną tabelę porównawczą statusów,
4. utrzymuje `strict_closure_claim_allowed=false` i `theory_closure_status=OPEN`.

## Honest boundary

`P1215` jest testem metodologicznym jakości bramkowania.
Nie stanowi dowodu rozładowania W4 i nie stanowi domknięcia teorii.
