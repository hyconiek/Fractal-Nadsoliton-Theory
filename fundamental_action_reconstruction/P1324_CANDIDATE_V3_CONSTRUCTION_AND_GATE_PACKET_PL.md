# P1324 — Candidate-v3 construction and gate packet (PL)

Status: `ROBUST_BUT_DEGENERATE_SELECTOR_CANDIDATE`
As of: `2026-05-12`
Depends on: `P1323`

## Cel
Po falsyfikacji `v2` zbudować `v3` i przeprowadzić podwójną bramkę:
1. robustność (brak kontrprzykładów),
2. informatywność selektora (nie-trivialne rozróżnianie znaków).

## Artefakt wykonawczy
- skrypt: `p1324_candidate_v3_construction_and_gate.py`
- raport: `generated/p1324_candidate_v3_construction_and_gate_report_v1.json`

## Wynik
- `negative_count = 0` (robustność na oknie testowym: tak),
- `sign_diversity = 1` (informatywność: nie),
- status: `ROBUST_BUT_DEGENERATE`.

## Decyzja profesorska
`S_sel_strict_v3` jest stabilny numerycznie, ale degeneruje się do
quasi-zawsze-dodatniego selektora. Taki obiekt nie jest dobrym kandydatem na
fizyczny selektor kierunku, bo nie niesie realnej mocy rozstrzygającej.

## Konsekwencja
- `v3` nie przechodzi bramki "admissible_as_strict_selector_candidate".
- `QW-2191` strict pozostaje `NOT_CLOSED`.
- Potrzebny `v4`: jednocześnie robust i informatywny.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że robustność sama wystarcza do selektora.

## Rekomendowany następny uczciwy krok
Uruchomić **P1325 balanced-v4 construction** z warunkiem brzegowym:
`negative_count = 0` **i** `sign_diversity = 2` na rozszerzonym oknie,
plus replay/adversarial jak w P1318.

## Dla laika
Zrobiliśmy kompas, który już się nie psuje na testach, ale ma inny problem:
prawie zawsze pokazuje to samo, więc niewiele rozstrzyga. Potrzebujemy wersji,
która jest jednocześnie stabilna i naprawdę odróżnia sytuacje.
