# P1389 Strict-only `ce6` Adaptive Boundary-Reweight Attempt (No Legacy Bridge) — PL

Status: `P1389_EXECUTED_ADAPTIVE_REWEIGHT_PARTIAL`
As of: `2026-05-13`

## Cel

Po `P1388` wykonać adaptacyjną re-wagę boundary patchy,
aby zredukować `sign_flip_rate` dla `ce6` poniżej progu
bez naruszania strict-only rygoru.

## Rygor

- `legacy_bridge_used = false`
- brak silent legacy transfer
- stałe kryterium: `sign_flip_rate <= epsilon_sign_v1`

## Metoda v1

1. Ustalić adaptacyjne wagi `w_boundary` zależne od lokalnej krzywizny overlap.
2. Wykonać rerun sygnatury znaku na tej samej rodzinie patchy.
3. Porównać `sign_flip_rate` z `P1388`.

## Wynik

`CE6_ADAPTIVE_REWEIGHT_STATUS := PARTIAL_IMPROVEMENT`

- `sign_flip_rate_before = 0.07` (z P1388)
- `sign_flip_rate_after = 0.052`
- poprawa jest, ale nadal `0.052 > epsilon_sign_v1=0.05`.

## Konsekwencja

- `ce6` formalnie nadal `UNRESOLVED`.
- `L_B1_03_EXPORT_STATUS` pozostaje `NOT_EXPORTED`.
- `B1` pozostaje `OPEN`.

## Decyzja profesorska

Następny krok: `P1390_STRICT_ONLY_CE6_EPSILON_ROBUSTNESS_SWEEP`
- sweep wokół granicy `epsilon_sign_v1`,
- test stabilności wyniku na perturbacjach siatki i wag,
- eksport: `robust_fail` albo `robust_pass` bez zmiany gate'ów.

## Omówienie dla laika

To jak precyzyjne strojenie instrumentu: jesteśmy już bardzo blisko tonu docelowego,
ale jeszcze minimalnie „fałszuje”, więc potrzebna jest ostatnia seria dokładnych testów.
