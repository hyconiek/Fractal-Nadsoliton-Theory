# P1409 — STRICT-ONLY PC5 EDGE ROBUSTNESS SWEEP (PL)

## Cel kroku
Po edge-pass z `P1408` testujemy odporność `PC5` na perturbacje siatki/wag
bez bridge do legacy i bez zmiany progów.

## Progi (zamrożone)
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`

## Wynik sweepu
- `num_trials = 20`
- `sign_pass_count = 20/20`
- `drift_pass_count = 14/20`
- `drift_min = 0.037`
- `drift_median = 0.040`
- `drift_max = 0.043`
- `PC5_EDGE_ROBUSTNESS_VERDICT := ROBUST_FAIL`

## Wniosek rygorystyczny
Mimo silnej poprawy jakości, brak pełnej stabilności robust-pass:
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1410_STRICT_ONLY_PC5_DRIFT_LOCAL_OBSTRUCTION_EXPORT_AND_PC6_TRIGGER` i przejść noncyklicznie do nowej klasy `PC6`.

## Omówienie dla laika
Nowa wersja działa wyraźnie lepiej, ale nadal nie przechodzi wszystkich wariantów testu odporności. To znaczy, że trzeba zrobić kolejny skok jakościowy, a nie tylko drobne strojenie.
