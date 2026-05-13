# P1387 Strict-only `ce6` Reparametrization Obstruction-or-Discharge (No Legacy Bridge) — PL

Status: `P1387_EXECUTED_CE6_CORE_ANALYSIS_OBSTRUCTION_EXPORTED`
As of: `2026-05-13`

## Cel

Po `P1386` rozstrzygnąć ostatni aktywny blocker (`ce6`) przez:
- theorem-level discharge,
  **albo**
- jawny eksport obstruction core,
w torze strict-only bez legacy bridge.

## Rygor

- `legacy_bridge_used = false`
- zero silent role transfer
- no-false-pass: brak sztucznego zamknięcia `L-B1-03`

## Analiza `ce6`

`ce6` dotyczy odporności na reparametryzację overlap przy stałym blocker-cut.

Wykonane testy v1:
1. perturbacje parametrów overlap (`rho`, `sigma`) w dopuszczalnym zakresie,
2. porównanie dryfu selector-compatibility,
3. test stabilności znaku operatora na granicy patchy.

## Wynik

`CE6_VERDICT := OBSTRUCTION_EXPORTED_V1`

- znaleziono niestabilność znaku na podobszarze granicznym,
- drift przekracza bound w części trajektorii,
- theorem-level discharge na dziś: `NOT_ACHIEVED`.

## Konsekwencja

- `L_B1_03_EXPORT_STATUS` pozostaje `NOT_EXPORTED`.
- `B1` pozostaje `OPEN`.
- Program pozostaje uczciwy: jawna przeszkoda zamiast pozornego PASS.

## Decyzja profesorska

Następny krok: `P1388_STRICT_ONLY_CE6_SIGN_STABILITY_PATCH_REFINEMENT_ATTEMPT`
- lokalna rafinacja patch-boundary,
- formalny sign-stability lemma candidate,
- rerun `ce6` z tym samym no-false-pass gate.

## Omówienie dla laika

To jak znalezienie mikropęknięcia w kluczowym elemencie konstrukcji.
Nie udajemy, że go nie ma — zapisujemy dokładnie gdzie jest,
i dopiero potem projektujemy poprawkę.
