# P1539 S489 TB_LEM4 Noise-Stability Theorem-Level Checkpoint Packet (No Legacy Bridge)

Status: `P1539_EXECUTED_TB_LEM4_NOISE_STABILITY_THEOREM_LEVEL_CHECKPOINT_PROVISIONAL`
As of: `2026-05-14`

## Cel

Następny uczciwy krok po `P1538`:

- formalnie sprawdzić `TB_LEM_4_NOISE_STABILITY`,
- ocenić czy noise-branch nie zmienia decyzji tie-break,
- utrzymać strict-only i brak legacy bridge.

## Zakres

`S489` uruchamia rodzinę perturbacji noise-branch i mierzy zgodność decyzji
z runem bazowym.

## Kontrakt wyjścia

- `base_decision`,
- `noise_scenarios_results`,
- `noise_stability_pass` (bool),
- `tb_lem4_status_update` in `{theorem_level_candidate, partial, open}`,
- `qw2191_closed=false`.

## PASS/FAIL

PASS jeśli wszystkie scenariusze noise zachowują tę samą decyzję tie-break.

FAIL jeśli choć jeden scenariusz zmienia decyzję.
