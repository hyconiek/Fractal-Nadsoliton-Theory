# P1767 / S717 — Sequencing: Bianchi/Ward -> BRST -> Cutkosky

Status: `P1767_S717_STRICT_GATE_SEQUENCING_CONTRACT_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Wprowadzono jawny kontrakt sekwencjonowania theorem-gates po aktualizacji
state-vector (`P1766`):

1. `G_BW`: Bianchi/Ward divergence gate,
2. `G_BRST`: BRST nilpotency gate,
3. `G_CUT`: Cutkosky unitarity gate.

Sekwencja jest twarda i strict-only: brak możliwości pominięcia `G_BW`.

## Co zostało dowiedzione

1. Ustalono formalną regułę promocji:
   **bez `G_BW: PASS_ZERO` nie wolno promować bramek QG theorem-level**.
2. Każda bramka ma zdefiniowane:
   - wymagane wejścia,
   - dopuszczalne wyjścia,
   - status blokady.

## Co nadal jest OPEN

1. `G_BW` jest OPEN (brak komponentowego divergence/residual execution).
2. `G_BRST` jest OPEN i zablokowane przez `G_BW`.
3. `G_CUT` jest OPEN i zablokowane przez `G_BW + G_BRST`.

## Ryzyka false-pass

1. Formalne „przeskoczenie” do BRST/Cutkosky bez zamknięcia Bianchi/Ward.
2. Mylenie gotowości kontraktowej z dowodem fizycznym.
3. Niejawne mieszanie klas `REDUCED/NONPROXY` lub `LOCAL/GLOBAL`.

## Następny uczciwy krok

Uruchomić `G_BW` na tej samej rodzinie teł:

- policzyć `EL_g-E_{μν}` komponentowo,
- policzyć `∇_μ(E_total^{μν})`,
- wydać wyłącznie `PASS_ZERO` albo `OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

Dopiero po `PASS_ZERO` wolno otworzyć formalne wykonanie `G_BRST`.

## Dla laika

To jest jak kolejka testów bezpieczeństwa: najpierw test podstawowej spójności
równań, potem test symetrii kwantowej, a na końcu test unitarności. Jeśli
pierwszy test nie przejdzie, następne nie mają sensu.
