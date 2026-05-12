# P1322 — Internal selector-source construction packet (PL)

Status: `CANDIDATE_INTERNAL_SELECTOR_SOURCE_V2_EXPORTED_NOT_DISCHARGED`
As of: `2026-05-12`
Depends on: `P1321`

## Cel
Wykonać rekomendowany krok z `P1321`: zbudować jawny kandydat
`S_sel_strict_v2` niezależny od `z2/eps` i sprawdzić, czy na obecnym zbiorze
C1–C4 daje stabilny znak kierunku.

## Artefakt wykonawczy
- skrypt: `p1322_internal_selector_source_construction_probe.py`
- raport: `generated/p1322_internal_selector_source_construction_probe_report_v1.json`

## Kandydat
`S_sel_strict_v2(phase, amp) = sign((phase-0.30) + 0.5*(amp-0.50))`

Własność konstrukcyjna: brak zależności od `z2/eps`.

## Wynik lokalny
- `uniqueness_on_current_dataset = true`
- `distinct_branch_signs = [+1]`
- status: `CANDIDATE_PASSES_LOCAL_DATASET`
- nadal: `strict_core_selector_source_exported = false`
- nadal: `qw2191_strict_status = NOT_CLOSED`

## Decyzja profesorska
To jest **konstrukcja kandydata**, nie domknięcie.

Mamy teraz konkretną hipotezę `S_sel_strict_v2`, która lokalnie eliminuje
zależność od residual slotu, ale brakuje jeszcze:
1. theorem-level derivation strict-side,
2. globalnego testu falsyfikacji poza C1–C4,
3. dowodu odporności na pełne klasy perturbacji.

## Konsekwencja
- Program przechodzi z etapu "brak kandydata" do etapu "kandydat jawny".
- `QW-2191` strict nadal otwarte do czasu discharge punktów 1–3.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że kandydat jest już prawem wyprowadzonym z teorii.

## Rekomendowany następny uczciwy krok
Uruchomić **P1323 falsification-and-robustness battery**:
przetestować `S_sel_strict_v2` na rozszerzonym zbiorze i perturbacjach,
a następnie sformalizować warunki, pod którymi można awansować go na
`strict_core_selector_source_exported = true`.

## Dla laika
To jak zbudowanie nowego kompasu, który na aktualnych testach pokazuje
stabilnie jeden kierunek i nie zależy od problematycznego pokrętła `z2/eps`.
Ale zanim uznamy go za "oficjalny", trzeba go sprawdzić w trudniejszych
warunkach i udowodnić, że działa nie tylko na małej próbce.
