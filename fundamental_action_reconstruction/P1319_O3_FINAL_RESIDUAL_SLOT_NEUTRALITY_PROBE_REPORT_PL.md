# P1319 — O3-final residual-slot neutrality probe report (strict-only, PL)

Status: `O3_FINAL_GATE_EXECUTED_NEUTRALITY_NOT_PROVEN`
As of: `2026-05-12`
Depends on: `P1318`

## Cel
Wykonać finałową bramkę O3: sprawdzić, czy residual slot `open(Z2/eps)` jest
neutralny względem wyboru kierunku (wtedy można domykać), czy nadal może
produkować alternatywną legalną gałąź (wtedy nie wolno domykać).

## Artefakt wykonawczy
- skrypt: `p1319_o3_final_residual_slot_neutrality_probe.py`
- raport: `generated/p1319_o3_final_residual_slot_neutrality_probe_report_v1.json`

## Wynik
- `distinct_branch_signs = [-1, +1]`
- `neutrality_holds = false`
- `verdict = NEUTRALITY_NOT_PROVEN`
- `nonuniqueness_residual_status = OPEN`
- `qw2191_strict_status = NOT_CLOSED`

## Decyzja profesorska
To jest rozstrzygnięcie negatywne dla próby domknięcia na obecnym stanie:
residual slot nadal może odwracać kierunek w legalnym obszarze parametrów,
więc unikalność strict-core nie została udowodniona.

## Konsekwencja naukowa
- Uczciwy status końcowy: `QW-2191` nadal otwarte w strict lane.
- Dalsze claimy closure bez nowego prawa selektora byłyby nienaukowe.

## Czego dokument nie twierdzi
- Nie twierdzi ToE closure.
- Nie twierdzi, że problem jest nierozwiązywalny.
- Twierdzi tylko, że na obecnym formalizmie domknięcie nie przeszło.

## Rekomendowany następny uczciwy krok
Zaprojektować i wyeksportować **nowy internal strict selector law**
(`S_sel_strict_v2`) eliminujący zależność znaku od `open(Z2/eps)` albo
udowodnić theorem-level constraint, który czyni gałąź `z2=-1` niedopuszczalną.

## Dla laika
Sprawdziliśmy "ostatnią śrubkę" i okazało się, że w pewnych legalnych
ustawieniach kompas może się odwrócić. To znaczy: jeszcze nie wolno ogłosić,
że teoria jest domknięta — trzeba najpierw zablokować tę możliwość.
