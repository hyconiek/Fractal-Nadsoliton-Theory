# P1405 — STRICT-ONLY PC4 EDGE ROBUSTNESS SWEEP (PL)

## Cel kroku
Po edge-pass z `P1404` wykonujemy sweep odporności `PC4` na perturbacje siatki/wag,
bez bridge do legacy i bez zmiany progów.

## Progi (zamrożone)
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`

## Wynik sweepu
- `num_trials = 20`
- `sign_pass_count = 20/20`
- `drift_pass_count = 9/20`
- `drift_min = 0.038`
- `drift_median = 0.041`
- `drift_max = 0.045`
- `PC4_EDGE_ROBUSTNESS_VERDICT := ROBUST_FAIL`

## Wniosek rygorystyczny
Edge-pass nie generalizuje stabilnie. Dlatego:
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1406_STRICT_ONLY_PC4_DRIFT_LOCAL_OBSTRUCTION_EXPORT_AND_PC5_TRIGGER`: formalnie wyeksportować `PC4-DRIFT-v1`, zamknąć pętlę `PC4` noncyklicznie i przejść do nowej klasy `PC5`.

## Omówienie dla laika
W jednym ustawieniu model zdał test, ale po wielu perturbacjach wynik często wracał ponad próg. To oznacza, że poprawa nie jest jeszcze stabilna i trzeba zmienić klasę podejścia.
