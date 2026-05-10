# P1155 Failure-Stage Diagnostic Refinement Packet

Status: `P1155_EXECUTED_FAILURE_STAGE_DIAGNOSTIC_REFINEMENT_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Po `P1154` (allow-failures) wykonujemy następny uczciwy krok:

nie tylko wiedzieć, że kandydat failuje, ale **na którym etapie** pipeline.

## Professor-level decision

Dodajemy diagnostykę failure-stage:

1. `P1151` exportuje `failed_step` (index + command),
2. `P1152` przenosi to do wyniku per-candidate jako `failure_stage`,
3. `P1153` pokazuje `failure_stage` w rankingu, aby priorytetyzować naprawy.

## Why this is strict-rigor

To redukuje niejednoznaczność metodologiczną: fail staje się lokalizowalny,
reprodukowalny i porównywalny między kandydatami.

## Honest boundary

`P1155` jest udoskonaleniem diagnostycznym procesu. Nie jest twierdzeniem
fizycznym i nie daje closure / `QW-2191` discharge.
