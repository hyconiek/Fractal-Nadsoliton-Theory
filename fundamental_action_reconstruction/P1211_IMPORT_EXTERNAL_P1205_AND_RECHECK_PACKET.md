# P1211 Import External P1205 And Recheck Packet

Status: `P1211_PREPARED_IMPORT_EXTERNAL_P1205_AND_RECHECK_NO_FALSE_PASS`
As of: `2026-05-11`

## Goal

Następny uczciwy krok: wpiąć rzeczywisty zewnętrzny artefakt CAS (`P1205`) do
repo i od razu wykonać kontrolowany re-check W4 bez ręcznych przepisań.

## Professor-level decision

Dodaję `p1211_import_external_p1205_and_recheck.py`, który:

1. importuje wskazany zewnętrzny JSON `P1205`,
2. aktualizuje lokalny `generated/p1205_...`,
3. liczy re-check statusu W4 zgodnie z polityką `P1210`.

## Honest boundary

`P1211` automatyzuje import i re-check, ale nie znosi guardrailu
`strict_closure_claim_allowed=false`.
