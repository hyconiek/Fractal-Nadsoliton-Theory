# P1400 — STRICT-ONLY PC3 EDGE ROBUSTNESS SWEEP (PL)

## Cel kroku
Po edge-pass z `P1399` wykonujemy test odporności (`robustness sweep`) dla klasy `PC3`:
- perturbacje siatki/wag,
- te same zamrożone progi,
- brak retuningu kryteriów,
- strict-only, bez bridge do legacy.

## Progi (zamrożone)
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`

## Wynik sweepu
- `num_trials = 20`
- `sign_pass_count = 20/20`
- `drift_pass_count = 7/20`
- `drift_min = 0.038`
- `drift_median = 0.041`
- `drift_max = 0.046`
- `PC3_EDGE_ROBUSTNESS_VERDICT := ROBUST_FAIL`

## Wniosek rygorystyczny
Mimo edge-pass w pojedynczym runie, wynik nie generalizuje stabilnie. Dlatego:
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1401_STRICT_ONLY_PC3_DRIFT_LOCAL_OBSTRUCTION_EXPORT_AND_PC4_TRIGGER`: formalnie wyeksportować lokalną przeszkodę `PC3-DRIFT-v1`, zamknąć pętlę `PC3` noncyklicznie i zaprojektować nową klasę `PC4`.

## Omówienie dla laika
W pojedynczym teście udało się zdać, ale kiedy sprawdziliśmy wiele wariantów — wynik często wracał ponad próg. To znaczy, że model jeszcze nie jest stabilny i trzeba uczciwie zmienić klasę podejścia, zamiast dalej stroić ten sam wariant.
