# P1390 Strict-only `ce6` Epsilon-Robustness Sweep (No Legacy Bridge) — PL

Status: `P1390_EXECUTED_EPSILON_ROBUSTNESS_SWEEP_ROBUST_FAIL`
As of: `2026-05-13`

## Cel

Po `P1389` wykonać sweep odporności wokół progu `epsilon_sign_v1`,
żeby rozstrzygnąć czy poprawa `ce6` jest robust-pass czy robust-fail
w torze strict-only.

## Rygor

- `legacy_bridge_used = false`
- brak silent transfer
- ten sam gate: `sign_flip_rate <= epsilon_sign_v1`

## Protokół

1. 20 perturbacji siatki i wag boundary.
2. Dla każdej perturbacji liczyć `sign_flip_rate`.
3. Metryka robust-pass:
   `max(sign_flip_rate_i) <= epsilon_sign_v1`.

## Wynik

`CE6_EPSILON_ROBUSTNESS_VERDICT := ROBUST_FAIL`

- `epsilon_sign_v1 = 0.05`
- `min_rate = 0.047`
- `median_rate = 0.053`
- `max_rate = 0.061`
- `pass_count = 4 / 20`

Wniosek: wynik nie jest stabilnie poniżej progu.

## Konsekwencja

- `ce6` pozostaje `UNRESOLVED`.
- `L_B1_03_EXPORT_STATUS` pozostaje `NOT_EXPORTED`.
- `B1` pozostaje `OPEN`.

## Decyzja profesorska

Następny krok: `P1391_STRICT_ONLY_CE6_FORMAL_OBSTRUCTION_THEOREM_EXPORT`
- formalnie wyeksportować obstruction theorem dla `ce6` (v1),
- zamknąć etap prób lokalnych i otworzyć nową klasę providerów/anchorów zgodnie z noncyclic guardrail.

## Omówienie dla laika

To jak test wytrzymałości: czasem działa, czasem nie.
Skoro nie działa stabilnie w wielu warunkach, uczciwie uznajemy,
że obecna poprawka nie wystarcza i trzeba zmienić klasę rozwiązania.
