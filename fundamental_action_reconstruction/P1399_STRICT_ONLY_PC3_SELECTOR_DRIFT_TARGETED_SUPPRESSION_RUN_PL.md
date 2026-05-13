# P1399 — STRICT-ONLY PC3 SELECTOR-DRIFT TARGETED SUPPRESSION RUN (PL)

## Cel kroku
Po `P1398` wykonujemy celowany run redukcji `selector_drift` w klasie `PC3`, przy zachowaniu:
- strict-only lane,
- bez bridge do legacy,
- bez zmiany zamrożonych progów `epsilon_sign_v1`, `epsilon_drift_v1`.

## Zakres strict-only
- `legacy_bridge_used = false`
- `inherits_from_pc2_loop = false`
- `threshold_retro_tuning = forbidden`

## Wejście
- `generated/p1398_strict_only_pc3_first_dual_metric_run_summary.json`

## Parametry progowe (bez zmian)
- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`

## Wynik runu celowanego
- `sign_flip_rate = 0.046` -> `PASS`
- `selector_drift = 0.040` -> `PASS` (edge-pass)
- `PC3_TARGETED_RUN_STATUS := PASS_EDGE`

## Ograniczenie rygoru
To jest pierwszy edge-pass, nie jeszcze robust-pass. Dlatego:
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1400_STRICT_ONLY_PC3_EDGE_ROBUSTNESS_SWEEP`: perturbacje siatki/wag wokół edge-pass i werdykt `ROBUST_PASS/ROBUST_FAIL` bez przesuwania progów.

## Omówienie dla laika
Udało się „dotknąć progu zdania” w obu testach, ale to dopiero pierwszy sukces. Teraz trzeba sprawdzić, czy wynik utrzyma się stabilnie w wielu wariantach, a nie tylko w jednym ustawieniu.
