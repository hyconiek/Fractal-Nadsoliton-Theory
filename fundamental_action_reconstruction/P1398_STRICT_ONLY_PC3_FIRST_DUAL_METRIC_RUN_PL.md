# P1398 — STRICT-ONLY PC3 FIRST DUAL-METRIC RUN (PL)

## Cel kroku
Wykonać pierwszy pełny run `PC3` na zamrożonych progach z `P1397`:
- `sign_flip_rate <= epsilon_sign_v1`
- `selector_drift <= epsilon_drift_v1`

bez retuningu kryteriów po fakcie i bez bridge do legacy.

## Zakres strict-only
- `legacy_bridge_used = false`
- `inherits_from_pc2_loop = false`
- brak transferu ról legacy -> strict

## Parametry testu (zamrożone)
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`

## Wynik pierwszego runu PC3
- `sign_flip_rate = 0.047` -> `PASS`
- `selector_drift = 0.043` -> `FAIL`
- `PC3_FIRST_RUN_STATUS := PARTIAL_PASS`

## Werdykt programu B1
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Znaczenie naukowe
Nowa klasa `PC3` poprawnie utrzymuje metrykę znaku pod progiem, ale nadal nie domyka metryki driftu selektora. To wskazuje na realny postęp względem pętli PC2, lecz nie daje prawa do theorem-level closure.

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1399_STRICT_ONLY_PC3_SELECTOR_DRIFT_TARGETED_SUPPRESSION_RUN`: celowana redukcja driftu przy utrzymaniu PASS dla `sign_flip_rate` i bez zmiany progów.

## Omówienie dla laika
Nowy model zdał pierwszy z dwóch kluczowych testów, ale drugi jeszcze nie. To znaczy: kierunek jest dobry, ale do pełnego „zaliczenia” brakuje jeszcze jednego stabilnego elementu.
