# P1403 — STRICT-ONLY PC4 FIRST DUAL-METRIC RUN (PL)

## Cel kroku
Wykonać pierwszy run `PC4` na zamrożonych progach i sprawdzić dwa warunki:
- `sign_flip_rate <= epsilon_sign_v1`
- `selector_drift <= epsilon_drift_v1`

bez bridge do legacy i bez retuningu progów.

## Progi
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`

## Wynik pierwszego runu PC4
- `sign_flip_rate = 0.045` -> `PASS`
- `selector_drift = 0.041` -> `FAIL`
- `PC4_FIRST_RUN_STATUS := PARTIAL_PASS`

## Wniosek rygorystyczny
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

Pierwszy run wskazuje poprawę stabilności znaku, ale nadal nie domyka metryki driftu selektora.

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1404_STRICT_ONLY_PC4_SELECTOR_DRIFT_TARGETED_SUPPRESSION_RUN`: celowana redukcja driftu przy utrzymaniu PASS dla metryki znaku i bez zmiany progów.

## Omówienie dla laika
Nowa wersja modelu zdała połowę egzaminu: jeden wskaźnik jest dobry, drugi jeszcze minimalnie za słaby. Teraz potrzebna jest precyzyjna poprawka tylko tam, gdzie nadal jest błąd.
