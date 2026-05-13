# P1394 Strict-only `PC2` Selector-Drift Suppression Run (No Legacy Bridge) — PL

Status: `P1394_EXECUTED_PC2_DRIFT_SUPPRESSION_PARTIAL`
As of: `2026-05-13`

## Cel

Po `P1393` wykonać celowany run redukcji `selector_drift`
w klasie `PC2`, przy zachowaniu `sign_flip_rate` poniżej progu
oraz bez zmiany gate’ów rygoru.

## Rygor

- `legacy_bridge_used = false`
- `provider_class = PC2`
- progi niezmienione: `epsilon_sign_v1=0.05`, `epsilon_drift_v1=0.04`
- brak nowych aksjomatów

## Metoda v1

1. Wzmocnienie tłumienia dryfu na overlapach (`drift_damping_lambda_v1`).
2. Utrzymanie anchor-map z `PC2` bez cofania do starej klasy providerów.
3. Jednoczesna kontrola metryk `sign_flip_rate` i `selector_drift`.

## Wynik

`PC2_DRIFT_SUPPRESSION_STATUS := PARTIAL_PASS`

- `sign_flip_rate = 0.048` -> **PASS**
- `selector_drift = 0.041` -> **FAIL** (bardzo blisko progu)

## Konsekwencja

- `L_B1_03_EXPORT_STATUS` pozostaje `NOT_EXPORTED`.
- `B1` pozostaje `OPEN`.
- aktywny blocker: niedomknięty `selector_drift`.

## Decyzja profesorska

Następny krok: `P1395_STRICT_ONLY_PC2_DRIFT_EPSILON_EDGE_ROBUSTNESS_RUN`
- test odporności przy brzegu `epsilon_drift_v1`,
- rozstrzygnięcie: czy niedobór 0.001 jest stabilny artefaktowo,
- jeśli fail robust, eksport lokalnej przeszkody PC2-drift-v1.

## Omówienie dla laika

To jak próba domknięcia drzwi, które prawie się zatrzaskują:
brakuje dosłownie milimetra. Trzeba sprawdzić, czy ten brak to przypadek,
czy stały problem konstrukcyjny.
