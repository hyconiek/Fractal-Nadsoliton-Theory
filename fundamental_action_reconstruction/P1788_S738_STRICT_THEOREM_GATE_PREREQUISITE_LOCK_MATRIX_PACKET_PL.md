# P1788 S738 Strict Theorem Gate Prerequisite Lock Matrix Packet (PL)

## Cel

Ustawić twardą macierz blokad dla bramek theorem-level (`BW`, `BRST`, `Cutkosky`),
aby żadna lokalna poprawa (`H1`, `EL_g`) nie była nadinterpretowana jako globalne domknięcie QG.

## Założenia strict-only

- Zero legacy bridges.
- Zero proxy shortcuts.
- Zero theorem-level PASS bez jawnych witnessów globalnych.

## Prerequisite lock matrix

1. `TG1_BIANCHI_WARD_GLOBAL`
   - wymaga: komponentowy divergence trace + zgodność na wspólnej rodzinie teł.
   - status teraz: `OPEN_LOCKED_BY_COMPONENTWISE_DIVERGENCE`.
2. `TG2_BRST_GLOBAL_NILPOTENCY`
   - wymaga: jawny nilpotency witness dla aktualnego nonproxy bundle.
   - status teraz: `OPEN_LOCKED_BY_NILPOTENCY_WITNESS`.
3. `TG3_CUTKOSKY_GLOBAL_UNITARITY`
   - wymaga: jawny cut-compatibility witness bez ghost-pole contradiction.
   - status teraz: `OPEN_LOCKED_BY_UNITARITY_WITNESS`.

## Dependency policy

- `H1` i `EL_g-E_{μν}` są **konieczne**, ale **niewystarczające** dla theorem-level closure.
- Lokalny PASS (`LOCAL/NONPROXY`) nie promuje automatycznie `GLOBAL/THEOREM`.
- Promotion jest dopuszczalny dopiero po jednoczesnym domknięciu `TG1+TG2+TG3`.

## Co zostało dowiedzione

- Repo ma już jawne sygnały `OPEN` dla BW/BRST/Cutkosky; brak podstaw do premature theorem PASS.
- Konieczna jest rozłączna walidacja globalnych witnessów dla każdej z 3 bramek.

## Ryzyka false-pass

1. Promotion po samym `H1`.
2. Promotion po samym lokalnym residual `EL_g-E_{μν}`.
3. Promotion bez jawnych globalnych witnessów BRST/Cutkosky.

## Następny uczciwy krok

Po wykonaniu unified `G1/G2/G5` uruchomić oddzielny global gate-pack:
`BW divergence witness -> BRST nilpotency witness -> Cutkosky unitarity witness`,
bez claimu PASS na etapach pośrednich.

## Objaśnienie dla laika

To trzy niezależne kontrole bezpieczeństwa. Nawet jeśli jedna wyjdzie dobrze,
cały system nie jest jeszcze domknięty, dopóki wszystkie trzy nie przejdą.
