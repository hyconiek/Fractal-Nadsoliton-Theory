# P1768 / S718 — Input-pack do wykonania G_BW komponentowo

Status: `P1768_S718_READY_FOR_COMPONENTWISE_EXECUTION_ATTEMPT_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Zamiast kolejnego scaffoldu ogólnego przygotowano **konkretny input-pack wykonawczy**
dla `G_BW` (Bianchi/Ward gate), spięty z wcześniejszymi checkpointami:

- normalizacja indeksów i znaków (`P1716`),
- baza residualu `B1/B2/B3/C1/C2` + format współczynników (`P1727`),
- sekwencja bramek `G_BW -> G_BRST -> G_CUT` (`P1767`).

## Co zostało dowiedzione

1. Wejście do pierwszego realnego wykonania `G_BW` jest formalnie kompletne:
   - rodzina teł,
   - lock konwencji indeksowej,
   - baza residualu,
   - sloty wektora współczynników,
   - wymagane operatory nonproxy (`E_A^μ`, `E_H`, `EL_g`).
2. Wyjścia są twardo ograniczone:
   `PASS_ZERO` lub `OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Co nadal jest OPEN

1. Samo wykonanie komponentowe `EL_g-E_{μν}` jeszcze nie zostało uruchomione.
2. Brak policzonego `∇_μ(E_total^{μν})` trace.
3. BRST i Cutkosky nadal pozostają zablokowane przez brak wyniku `G_BW`.

## Ryzyka false-pass

1. Uznanie „input-pack gotowy” za dowód przejścia gate.
2. Brak pełnego wektora współczynników przy ogłaszaniu werdyktu.
3. Rozjechanie konwencji indeksowej między residualem i divergence-check.

## Następny uczciwy krok

Wykonać `G_BW` na tym input-packu i opublikować wynik z pełnym trace:

- `componentwise_EL_g_minus_E_munu_B1_B2_B3_C1_C2`,
- `componentwise_nabla_mu_E_total_munu_trace`,
- finalny werdykt: `PASS_ZERO` albo `OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Dla laika

To jak przygotowanie kompletnego zestawu laboratoryjnego przed eksperymentem:
wszystkie narzędzia są gotowe i zgodne, ale wynik naukowy pojawi się dopiero po
realnym wykonaniu pomiaru.
