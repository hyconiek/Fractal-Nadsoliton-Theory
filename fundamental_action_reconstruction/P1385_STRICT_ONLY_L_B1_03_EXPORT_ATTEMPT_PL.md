# P1385 Strict-only `L-B1-03` Export Attempt (No Legacy Bridge) — PL

Status: `P1385_EXECUTED_LB103_EXPORT_ATTEMPT_WITH_COUNTEREXAMPLE_SWEEP`
As of: `2026-05-13`

## Cel

Po lokalnym PASS z `P1384` wykonać próbę formalnego eksportu
`L-B1-03` (selector compatibility discharge) w torze strict-only:
`F_Nadsoliton => L_SM + L_GR` bez mostu do legacy.

## Rygor

- `legacy_bridge_used = false`
- brak transferu ról legacy->strict
- no-false-pass: brak theorem export bez przejścia pełnego sweepu kontrprzykładów

## Protokół

1. Zdefiniować klasę testową kontrprzykładów `CE_v1 = {ce1..ce6}`.
2. Sprawdzić, czy lokalny selector-pass utrzymuje się bez ukrytego selector-closure.
3. Ocenić warunek eksportu `L-B1-03`:
   - `all_counterexamples_resolved = true`
   - `global_selector_compatibility_bound_verified = true`

## Wynik

`L_B1_03_EXPORT_STATUS := PARTIAL_NOT_EXPORTED`

- `resolved_counterexamples = 4/6`
- `unresolved = {ce4, ce6}`
- `global_selector_compatibility_bound_verified = false`

Wniosek: lokalny postęp utrzymany, ale brak prawa do theorem-level export.

## Konsekwencja

`G3` przechodzi z "OPEN_QW2191_GATED" do "PARTIAL_QW2191_NARROWED",
ale `B1` nadal pozostaje `OPEN`.

## Decyzja profesorska

Następny krok: `P1386_STRICT_ONLY_COUNTEREXAMPLE_CE4_CE6_RESOLUTION_ATTEMPT`
- osobny run dla `ce4` i `ce6`,
- jawny downgrade do obstruction jeśli którykolwiek pozostanie nierozwiązany.

## Omówienie dla laika

To jak egzamin z 6 zadań: 4 zrobione dobrze, 2 nadal nie.
To jest postęp, ale nie można jeszcze ogłosić pełnego zaliczenia.
