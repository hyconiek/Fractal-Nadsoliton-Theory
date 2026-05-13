# P1393 Strict-only `PC2` First Sign-Stability and Selector-Drift Run (No Legacy Bridge) — PL

Status: `P1393_EXECUTED_PC2_BASELINE_FIRST_RUN_PARTIAL_PASS`
As of: `2026-05-13`

## Cel

Wykonać pierwszy run na nowej klasie providerów `PC2` (z `P1392`) i sprawdzić,
czy nowy anchor poprawia jednocześnie:
1. `sign_flip_rate`,
2. `selector_drift`,
bez powrotu do legacy bridge.

## Rygor

- `legacy_bridge_used = false`
- `inherits_from_ce6_v1_loop = false`
- te same progi formalne: `sign_flip_rate <= epsilon_sign_v1`, `selector_drift <= epsilon_drift_v1`

## Parametry runu v1

- `epsilon_sign_v1 = 0.05`
- `epsilon_drift_v1 = 0.04`
- `provider_class_id = PC2_strict_boundary_anchor_v1`

## Wynik

`PC2_FIRST_RUN_STATUS := PARTIAL_PASS`

- `sign_flip_rate = 0.049` -> **PASS**
- `selector_drift = 0.046` -> **FAIL**

## Wniosek

Nowa klasa `PC2` usuwa blocker znaku (na tym runie),
ale nie domyka jeszcze driftu selektora.

## Konsekwencja

- `L_B1_03_EXPORT_STATUS` pozostaje `NOT_EXPORTED`.
- `B1` pozostaje `OPEN`.
- stara pętla `ce6-v1` pozostaje zamknięta; pracujemy wyłącznie w gałęzi PC2.

## Decyzja profesorska

Następny krok: `P1394_STRICT_ONLY_PC2_SELECTOR_DRIFT_SUPPRESSION_RUN`
- celowana redukcja `selector_drift` przy zachowaniu `sign_flip_rate` poniżej progu,
- brak zmiany progów i brak aksjomatów ad hoc.

## Omówienie dla laika

Nowa metoda rozwiązała pierwszy problem, ale drugi jeszcze nie.
To jak naprawa auta: silnik już pracuje dobrze, ale trzeba jeszcze wyregulować układ kierowniczy.
