# P1770 / S720 — Sanitization notacji tensorowej i kontrakt retry dla G_BW

Status: `P1770_S720_NOTATION_SANITIZED_GBW_RETRY_CONTRACT_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Po `P1769` wykryto praktyczny blocker techniczny: niejednoznaczny zapis
indeksów krzywizny (ryzyko błędnej interpretacji przez backend/parser).

W `P1770` dodano:

1. oczyszczony, parser-safe lock notacji tensorowej (`ASCII-safe`),
2. kontrakt retry dla `G_BW` bez zmiany polityki werdyktu,
3. jawne utrzymanie tej samej rodziny teł i bazy `B1/B2/B3/C1/C2`.

## Co zostało dowiedzione

1. Część obstrukcji z `P1769` miała komponent notacyjno-implementacyjny,
   a nie tylko czysto fizyczny.
2. Można uczciwie powtórzyć `G_BW` na tych samych danych po naprawie notacji,
   bez naruszania strict-only route.

## Co nadal jest OPEN

1. Finalny werdykt `G_BW` po retry (`PASS_ZERO` lub `OBSTRUCTION_WITH_DIVERGENCE_TRACE`).
2. Odblokowanie `G_BRST` i `G_CUT`.
3. Theorem-level closure (renormalization/BRST/Cutkosky/background-independence).

## Ryzyka false-pass

1. Zmiana notacji nie może być traktowana jako dowód fizycznego domknięcia.
2. Retry nie może zmieniać kryteriów PASS/OBSTRUCTION.
3. Nie wolno zmieniać rodziny teł ani bazy residualu „po cichu”.

## Następny uczciwy krok

Uruchomić retry `G_BW` z oczyszczoną notacją i opublikować pełny trace:

- `componentwise_EL_g_minus_E_munu_B1_B2_B3_C1_C2`,
- `componentwise_nabla_mu_E_total_munu_trace`,
- finalny werdykt: `PASS_ZERO` albo `OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Dla laika

To jak poprawa zapisu wzoru przed włożeniem go do kalkulatora: sama poprawka
zapisu nie daje wyniku, ale usuwa błąd techniczny, który mógł blokować
poprawne liczenie.
