# P1404 — STRICT-ONLY PC4 SELECTOR-DRIFT TARGETED SUPPRESSION RUN (PL)

## Cel kroku
Po `P1403` wykonujemy celowaną korektę tylko dla metryki `selector_drift` w klasie `PC4`,
bez zmiany progów i bez bridge do legacy.

## Progi (zamrożone)
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`

## Wynik runu celowanego
- `sign_flip_rate = 0.044` -> `PASS`
- `selector_drift = 0.040` -> `PASS` (edge-pass)
- `PC4_TARGETED_RUN_STATUS := PASS_EDGE`

## Ograniczenie rygoru
To jest pierwszy edge-pass i jeszcze nie robust-pass. Zatem:
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1405_STRICT_ONLY_PC4_EDGE_ROBUSTNESS_SWEEP`: test odporności na perturbacje siatki/wag i jawny werdykt `ROBUST_PASS/ROBUST_FAIL` bez przesuwania progów.

## Omówienie dla laika
Udało się poprawić dokładnie ten element, który wcześniej zawodził. Teraz trzeba sprawdzić, czy to stabilna poprawa, czy tylko jednorazowy dobry wynik.
