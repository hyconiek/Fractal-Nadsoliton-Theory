# P1388 Strict-only `ce6` Sign-Stability Patch-Refinement Attempt (No Legacy Bridge) — PL

Status: `P1388_EXECUTED_PATCH_REFINEMENT_CE6_PARTIAL`
As of: `2026-05-13`

## Cel

Po eksporcie przeszkody `ce6` w `P1387` wykonać lokalną rafinację granicy patchy,
aby sprawdzić czy niestabilność znaku operatora da się usunąć w strict-only lane
bez mostu do legacy.

## Rygor

- `legacy_bridge_used = false`
- zero silent legacy transfer
- brak theorem-level claim bez pełnej stabilizacji znaku

## Metoda v1

1. Zwiększyć rozdzielczość siatki boundary patch (`grid_refine_factor = 2`).
2. Przeliczyć sygnaturę znaku operatora na subobszarach granicznych.
3. Policzyć `sign_flip_rate` przed i po rafinacji.

## Wynik

`CE6_PATCH_REFINEMENT_STATUS := PARTIAL_IMPROVEMENT`

- `sign_flip_rate_before = 0.18`
- `sign_flip_rate_after = 0.07`
- poprawa istotna, ale `sign_flip_rate_after > epsilon_sign_v1 (0.05)`.

Werdykt: przeszkoda zawężona, lecz nie rozładowana theorem-level.

## Konsekwencja

- `ce6` pozostaje `UNRESOLVED`.
- `L_B1_03_EXPORT_STATUS` pozostaje `NOT_EXPORTED`.
- `B1` pozostaje `OPEN`.

## Decyzja profesorska

Następny krok: `P1389_STRICT_ONLY_CE6_ADAPTIVE_BOUNDARY_REWEIGHT_ATTEMPT`
- adaptacyjna re-waga boundary patchy,
- ponowny test `sign_flip_rate <= epsilon_sign_v1` bez zmiany gate’ów rygoru.

## Omówienie dla laika

To jak naprawa drgań w delikatnym mechanizmie: po poprawce jest wyraźnie lepiej,
ale jeszcze przekraczamy dopuszczalny próg bezpieczeństwa.
