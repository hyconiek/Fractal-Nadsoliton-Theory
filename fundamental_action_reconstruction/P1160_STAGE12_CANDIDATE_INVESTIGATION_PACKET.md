# P1160 Stage1/2 Candidate Investigation Packet

Status: `P1160_EXECUTED_STAGE12_CANDIDATE_INVESTIGATION_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Na prośbę "badaj kandydatów" wykonuję następny uczciwy krok:

1. utworzyć kontrolne kandydaty z awarią na stage1 i stage2,
2. uruchomić registry + explorer,
3. sprawdzić czy `P1159` faktycznie działa na realnych targetach.

## Professor-level decision

Dodano jawny mechanizm testowy `force_fail_stage` w `P1151` (tylko do
kontrolowanych eksperymentów metodologicznych), aby tworzyć replikowalne
kandydaty stage1/stage2 i badać działanie narzędzi naprawczych.

## Added candidate artifacts

- `generated/p1160_candidate_stage1_fail.json`
- `generated/p1160_candidate_stage2_fail.json`
- `generated/p1160_candidate_registry_stage12.json`

## Result

1. `P1152` na nowym rejestrze: `total=2, passed=0, failed=2`,
2. `P1159` wykrył targety stage1/2 (`target_count=2`) i wykonał 2 próby,
3. raport jest wygenerowany jawnie w
   `generated/p1159_stage12_repair_explorer_summary.json`.

## Honest boundary

To nadal jest badanie metodologiczne kandydatów i narzędzi naprawy.
Brak claimu closure i brak discharge `QW-2191`.
