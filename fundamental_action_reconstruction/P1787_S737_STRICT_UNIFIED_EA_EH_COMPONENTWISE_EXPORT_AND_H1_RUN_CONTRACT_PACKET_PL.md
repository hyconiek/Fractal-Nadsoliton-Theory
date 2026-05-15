# P1787 S737 Strict Unified EA/EH Componentwise Export And H1 Run Contract Packet (PL)

## Cel

Zdefiniować jeden spójny strict-only run-contract, który łączy:

1. komponentowy export `E_A^μ` (nonproxy),
2. komponentowy export `E_H` (nonproxy),
3. test `H1`: `δE_A^μ/δH - δE_H/δA_μ`,

na **tej samej rodzinie teł** i w **tej samej konwencji indeksowej**.

## Dlaczego to jest krok o najwyższej wartości

Aktualny bottleneck nie wynika z braku kolejnych scaffoldów,
lecz z braku *unified execution contract* dla `G1+G2+G5`.
Bez tego łatwo o false-pass przez niespójne tła lub konwencje indeksowe.

## Strict execution contract

- `route_policy`: strict-only, zero legacy bridge.
- `artifact_policy`: nonproxy-only dla eksportów wejściowych.
- `background_policy`: jedna jawnie zamrożona rodzina teł dla wszystkich 3 kroków.
- `index_policy`: jedna jawnie zamrożona konwencja indeksowa.
- `pass_policy`: PASS tylko przy jawnym residual vector = 0 i witness ledger.

## Gate mapping

- `G1_EA_NONPROXY_EXPLICIT_EXPORT` -> wymagany output: komponentowy ledger `E_A^μ`.
- `G2_EH_NONPROXY_EXPLICIT_EXPORT` -> wymagany output: komponentowy ledger `E_H`.
- `G5_H1_4D_WEAK_FORM` -> wymagany output: komponentowy ledger różnicy
  `δE_A^μ/δH - δE_H/δA_μ`.

## Status discipline

- Jeśli brak pełnej ekspansji komponentowej: `OPEN_COMPONENTWISE_REQUIRED`.
- Jeśli residual policzony, ale niezerowy: `OPEN_OBSTRUCTION_WITH_TRACE`.
- Tylko przy pełnym zerze i kompletnym witness: `PASS_STRICT_LOCAL_H1_ZERO`.
- Niezależnie od wyniku lokalnego: theorem-level QG pozostaje `OPEN`
  do czasu BRST/Cutkosky/Bianchi/Ward domknięcia.

## Ryzyka false-pass

1. Mieszanie różnych rodzin teł między `E_A^μ` i `E_H`.
2. Mieszanie różnych konwencji indeksowych między eksportem i H1.
3. Deklaracja PASS bez publikacji residual ledgers.

## Następny uczciwy krok

Uruchomić pojedynczy generator-checkpoint, który scala statusy `P1764/P1760/P1761/P1762/P1786`
i wymusza jeden run-order `G1 -> G2 -> G5` bez przeskoków.

## Objaśnienie dla laika

To jak test laboratoryjny: oba pomiary i porównanie muszą być wykonane
na tym samym sprzęcie i w tych samych warunkach, inaczej wynik nie jest wiarygodny.
