# P1395 Strict-only `PC2` Drift Epsilon-Edge Robustness Run (No Legacy Bridge) — PL

Status: `P1395_EXECUTED_PC2_DRIFT_EDGE_SWEEP_ROBUST_FAIL`
As of: `2026-05-13`

## Cel

Po `P1394` rozstrzygnąć, czy niedobór `selector_drift` o ~0.001 jest stabilny,
czy artefaktowy, przez test odporności przy brzegu `epsilon_drift_v1`.

## Rygor

- `legacy_bridge_used = false`
- progi bez zmian (`epsilon_sign_v1=0.05`, `epsilon_drift_v1=0.04`)
- brak aksjomatów ad hoc

## Protokół

1. 24 perturbacje siatki i parametrów tłumienia dryfu.
2. Dla każdej perturbacji liczyć `selector_drift_i` i `sign_flip_rate_i`.
3. Warunek robust-pass drift:
   `max(selector_drift_i) <= epsilon_drift_v1`.

## Wynik

`PC2_DRIFT_EDGE_VERDICT := ROBUST_FAIL`

- `drift_min = 0.038`
- `drift_median = 0.042`
- `drift_max = 0.047`
- `drift_pass_count = 5 / 24`
- `sign_pass_count = 24 / 24`

Wniosek: problem driftu jest strukturalny w tej konfiguracji PC2-v1,
nie tylko losowy artefakt pojedynczego runu.

## Konsekwencja

- `L_B1_03_EXPORT_STATUS` pozostaje `NOT_EXPORTED`.
- `B1` pozostaje `OPEN`.
- eksport lokalnej przeszkody `PC2-DRIFT-v1` uzasadniony.

## Decyzja profesorska

Następny krok: `P1396_STRICT_ONLY_PC2_DRIFT_LOCAL_OBSTRUCTION_EXPORT_AND_PC3_DESIGN_TRIGGER`
- formalnie wyeksportować `PC2-DRIFT-v1`,
- uruchomić projekt nowej klasy `PC3` (noncyclic anchor class).

## Omówienie dla laika

To jak test samochodu na różnych drogach: na części tras działa dobrze,
ale na wielu nadal wypada poza normę. To znaczy, że trzeba zmienić konstrukcję,
a nie tylko dalej „dostrajać” to samo ustawienie.
