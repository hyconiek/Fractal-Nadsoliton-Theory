# P1159 Stage1/2 Repair Explorer Packet

Status: `P1159_EXECUTED_STAGE12_REPAIR_EXPLORER_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać następny uczciwy krok po `P1158`: rozszerzyć pętlę napraw na kandydatów
odpadających na etapach probe (`failure_stage_index in {1,2}`).

## Professor-level decision

Wprowadzam eksplorator `P1159`, który:

1. czyta wyniki rejestru (`P1152`),
2. wybiera kandydatów z fail stage 1/2,
3. tworzy konserwatywne warianty (tylko dla jawnych pól hint),
4. rerun `P1151` i zapisuje wynik.

Jeśli nie ma targetów stage1/2, raportuje jawne `skipped_reason` zamiast
udawania postępu.

## Artifact

- script:
  `p1159_stage12_repair_explorer.py`
- summary:
  `generated/p1159_stage12_repair_explorer_summary.json`

## Current result

Na aktualnym rejestrze:

```text
target_count = 0
attempt_count = 0
skipped_reason = no stage1/stage2 failures in current registry
```

To jest poprawny wynik rygorystyczny: brak sztucznego „naprawiania” bez danych.

## Honest boundary

`P1159` nie domyka teorii i nie rozwiązuje `QW-2191`.
