# P1183 Candidate Promotion Gate Packet

Status: `P1183_EXECUTED_CANDIDATE_PROMOTION_GATE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Zamienić decyzję o promocji kandydata na jawny, automatyczny gate oparty o
wcześniejsze artefakty rygoru.

## Professor-level decision

Dodaję `p1183_candidate_promotion_gate.py`, który promuje kandydata tylko gdy:

1. `P1182 band_lock_pass == true`,
2. `P1152` ma `require_safe_region == true`,
3. `P1152` ma poprawne shortlist consistency,
4. `P1152` ma `failed == 0`.

## Honest boundary

`P1183` nie daje closure i nie rozwiązuje `QW-2191`; to automatyzacja decyzji
operacyjnej.
