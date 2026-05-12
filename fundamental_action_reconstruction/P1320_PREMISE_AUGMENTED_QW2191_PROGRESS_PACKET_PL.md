# P1320 — Premise-augmented QW-2191 progress packet (PL)

Status: `PREMISE_AUGMENTED_CLOSURE_OBTAINED_NON_STRICT`
As of: `2026-05-12`
Depends on: `P1319`

## Cel
Podjąć kolejny uczciwy krok po negatywnym wyniku strict-final (`P1319`):
sprawdzić, czy jawna premise selektora może domknąć `QW-2191` w trybie
**non-strict** (zgodnie z guardrail: jawnie oznaczone, bez udawania strict closure).

## Premise
`P_sel_v1`: gałąź selektora jest dopuszczalna tylko gdy `z2 = +1`.

To jest jawne symmetry-breaking premise, a nie twierdzenie strict-core.

## Artefakt wykonawczy
- skrypt: `p1320_premise_augmented_selector_law_probe.py`
- raport: `generated/p1320_premise_augmented_selector_law_probe_report_v1.json`

## Wynik
- pod premise `P_sel_v1` wszystkie dopuszczalne próbki mają jeden znak kierunku,
- `uniqueness_under_premise = true`,
- `qw2191_status_under_premise = CLOSED_NON_STRICT`,
- jednocześnie: `qw2191_strict_status = NOT_CLOSED`.

## Decyzja profesorska
1. Osiągnięto realny postęp metodologiczny: closure warunkowe jest możliwe,
   jeśli jawnie przyjmujemy premise selektora.
2. Nie osiągnięto strict-core closure.
3. Każde raportowanie musi utrzymać etykietę: `NON_STRICT_PREMISE_AUGMENTED`.

## Konsekwencja
- Można kontynuować rozwój przewidywań/falsifikacji pod `P_sel_v1`.
- Nie można jeszcze twierdzić, że `QW-2191` jest domknięte bezwarunkowo strict.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie ukrywa premise jako "wynik automatyczny".

## Rekomendowany następny uczciwy krok
Rozpocząć **P1321 strictification program**: spróbować wyprowadzić `P_sel_v1`
wewnętrznie ze strict-side źródła (bez zewnętrznej premise). Jeśli się nie da,
utrzymać status trwały: `CLOSED_NON_STRICT / OPEN_STRICT`.

## Dla laika
To jak dodanie jednoznacznej reguły "ten typ kompasu jest legalny, tamten nie".
Po tej regule kierunek robi się jednoznaczny. Ale dopóki ta reguła nie wynika
sama z teorii podstawowej, to jest to domknięcie warunkowe, a nie pełne.
