# P1766 / S716 — Aktualizacja state-vector + kontrakt bramki Bianchi/Ward

Status: `P1766_S716_STATE_VECTOR_UPDATED_WITH_BIANCHI_WARD_GATE_CONTRACT_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Po `P1764` (jawne `E_A^μ`, `E_H`) i `P1765` (jawne `EL_g^{μν}`) wykonano
aktualizację centralnego state-vectora forward/reverse.

Nowy stan:

- forward:
  - `K_strict`: locked,
  - `coefficient_map`: exported,
  - `full non-skeleton L_total`: exported,
  - `E_A^μ`, `E_H`, `EL_g^{μν}`: explicit nonproxy operator exports.
- reverse:
  - `H1`: open (componentwise required),
  - `EL_g-E_{μν}`: open (componentwise required),
  - `Bianchi/Ward`: open (wymaga divergence trace),
  - global Helmholtz: open.

## Co zostało dowiedzione

1. Pipeline strict-only ma już jawne operatory dla trzech kluczowych bramek:
   `E_A^μ`, `E_H`, `EL_g^{μν}`.
2. Można formalnie uruchomić kontrakt Bianchi/Ward bez zmiany klasy tła,
   bo zależności wejściowe zostały jawnie zdefiniowane.

## Co nadal jest OPEN

1. Brak komponentowego wyliczenia `EL_g-E_{μν}` na `B1/B2/B3/C1/C2`.
2. Brak wykonania divergence check: `∇_μ(E_total^{μν})=0`.
3. Brak theorem-level closure dla BRST/Cutkosky/renormalization/background-independence.

## Ryzyka false-pass

1. Uznanie state-vector update za dowód fizyczny (to tylko synchronizacja stanu).
2. Uznanie operatorowego eksportu za wynik residualny.
3. Ominięcie trace dla Bianchi/Ward i wydanie nieuzasadnionego PASS.

## Następny uczciwy krok

Wykonać dwa obliczenia na tym samym tle i tej samej konwencji indeksowej:

1. `EL_g-E_{μν}` komponentowo (`B1/B2/B3/C1/C2`),
2. `∇_μ(E_total^{μν})` jako bramka Bianchi/Ward.

Dopuszczalne wyjścia: tylko `PASS_ZERO` albo `OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Dla laika

To jak aktualizacja mapy projektu po zbudowaniu trzech kluczowych modułów.
Mapa pokazuje, że fundament jest gotowy, ale najważniejsze testy poprawności
muszą być teraz policzone liczbowo/symbolicznie krok po kroku.
