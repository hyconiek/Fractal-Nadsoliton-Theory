# P1778 — S728
## STRICT COMPONENTWISE H1 + METRIC JOINT EXECUTION LOCK PACKET (PL)

## Cel

Ustawić jeden, wspólny rygor wykonania dla dwóch krytycznych testów reverse-chain:
`H1` oraz `EL_g-E_{μν}` — na tej samej rodzinie teł, tej samej konwencji indeksowej
i tej samej bazie residualowej (`B1/B2/B3/C1/C2`).

## Technical progress

- Wyeksportowano formalny „joint execution lock” dla `R1 phase_1 + phase_2`.
- Jawnie zapisano zabronione ścieżki (mixed background, index drift, PASS bez residual witness).
- Utrzymano bramkę wejścia: brak startu dopóki `W1` nie uzyska `FULL_EXPORT`.

## Co zostało dowiedzione

1. Kontrakt spójnego uruchomienia H1+metric jest dostępny i audytowalny.
2. Aktualnie joint run jest poprawnie zablokowany przez niedomknięte W1.

## Co nadal jest OPEN

1. W1 FULL_EXPORT.
2. Componentwise witness H1.
3. Componentwise witness `EL_g-E_{μν}`.
4. Theorem-level Bianchi/Ward + BRST/Cutkosky.

## Ryzyka false-pass

- Uruchomienie H1 i metric na różnych tłach/konwencjach, a potem scalanie wyników jakby były porównywalne.
- Wydanie PASS bez jawnego residual-zero albo obstruction trace.

## Następny uczciwy krok

Domknąć W1 (`B2/B3` + normalizacja dywergencji `H(R2)`), następnie wykonać joint
H1+metric run zgodnie z lock-contract i opublikować wyłącznie:
`PASS_ZERO` albo `OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Wyjaśnienie dla laika

To jak ustawienie jednej procedury testowej dla dwóch krytycznych pomiarów,
żeby nie porównywać „jabłek z gruszkami”. Najpierw trzeba odblokować wejście,
a dopiero potem uruchomić oba testy razem i uczciwie pokazać wynik.
