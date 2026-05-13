# P1386 Strict-only Counterexample `ce4/ce6` Resolution Attempt (No Legacy Bridge) — PL

Status: `P1386_EXECUTED_TARGETED_CE4_CE6_RUN_PARTIAL`
As of: `2026-05-13`

## Cel

Po `P1385` wykonać celowany run rozwiązania nierozładowanych kontrprzykładów
`ce4` i `ce6` w torze strict-only dla ścieżki
`F_Nadsoliton => L_SM + L_GR`.

## Rygor

- `legacy_bridge_used = false`
- zakaz legacy-role transfer
- jawny downgrade do obstruction, jeśli którykolwiek CE pozostaje nierozwiązany

## Metoda

- `ce4`: test stabilności selector-compatibility przy rozszerzeniu nośnika.
- `ce6`: test odporności na reparametryzację overlap w tym samym blocker-cut.

## Wynik

`CE_RESOLUTION_STATUS := PARTIAL`

- `ce4 = resolved`
- `ce6 = unresolved`

`L_B1_03_EXPORT_STATUS` pozostaje `NOT_EXPORTED`.

## Konsekwencja

QW-2191 zawężony dalej, lecz nie rozładowany.
`B1` pozostaje `OPEN`.

## Decyzja profesorska

Następny krok: `P1387_STRICT_ONLY_CE6_REPARAMETRIZATION_OBSTRUCTION_OR_DISCHARGE`
- wyizolować minimalny obstruction core dla `ce6`,
- albo zamknąć go theorem-level i dopiero wrócić do exportu `L-B1-03`.

## Omówienie dla laika

Mieliśmy 2 nierozwiązane zadania. Jedno udało się zamknąć,
drugie nadal blokuje pełne zaliczenie. To dalej postęp,
ale uczciwie: finał jeszcze nie jest osiągnięty.
