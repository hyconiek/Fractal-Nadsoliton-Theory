# P1156 Repair-Priority Scoring Refinement Packet

Status: `P1156_EXECUTED_REPAIR_PRIORITY_SCORING_REFINEMENT_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Następny uczciwy krok po `P1155`: ranking ma wskazywać nie tylko jakość,
ale też priorytet naprawy kandydatów zablokowanych.

## Professor-level decision

Rozszerzam `P1153` o dwa jawne score:

1. `quality_score` (jak wcześniej: gate/pipeline/audit),
2. `repair_priority_score` dla faili:
   - im później kandydat odpada (`failure_stage.index` większy),
     tym wyższy priorytet naprawy.

To pozwala odróżnić „fundamentalnie słaby” fail od „prawie przechodzącego” faila.

## Outcome

Zaktualizowany summary `P1153` zawiera oba score i sortowanie:

1. najpierw po `quality_score`,
2. potem po `repair_priority_score`.

## Honest boundary

To nadal narzędzie metodologiczne. Brak claimu closure i brak discharge
`QW-2191`.
