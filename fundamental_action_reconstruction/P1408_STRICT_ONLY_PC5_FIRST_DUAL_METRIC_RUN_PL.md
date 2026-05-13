# P1408 — STRICT-ONLY PC5 FIRST DUAL-METRIC RUN (PL)

## Cel kroku
Wykonać pierwszy run `PC5` na zamrożonych progach:
- `sign_flip_rate <= epsilon_sign_v1`
- `selector_drift <= epsilon_drift_v1`

bez bridge do legacy i bez retuningu.

## Progi
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`

## Wynik pierwszego runu PC5
- `sign_flip_rate = 0.043` -> `PASS`
- `selector_drift = 0.040` -> `PASS` (edge-pass)
- `PC5_FIRST_RUN_STATUS := PASS_EDGE`

## Ograniczenie rygoru
To pierwszy edge-pass, jeszcze bez robust-pass:
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1409_STRICT_ONLY_PC5_EDGE_ROBUSTNESS_SWEEP` i sprawdzić stabilność na perturbacjach z werdyktem `ROBUST_PASS/ROBUST_FAIL`.

## Omówienie dla laika
To dobry sygnał: nowy model zdał oba testy w pierwszym podejściu. Ale naukowo trzeba jeszcze sprawdzić, czy wynik utrzyma się stabilnie w wielu wariantach.
