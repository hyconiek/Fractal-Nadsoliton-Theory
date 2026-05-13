# P1518 — S4.68 Strict Theorem Witness Compression And Channel Obligation Map Packet (PL)

Status: `P1518_EXECUTED_WITNESS_COMPRESSION_AND_CHANNEL_MAP`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok po `P1517`: skompresować `locked witness set`
do minimalnej bazy i jawnie zmapować ją na obowiązki kanałów
`F->LSM` oraz `F->LGR`.

## Decyzja profesorska

Nie zwiększamy zakresu twierdzenia. Redukujemy nośniki dowodu do minimum,
które nadal pokrywa oba kanały fizyczne strict-side.

## Wynik P1518

Publikujemy:

1. minimalną bazę witnessów,
2. mapę obowiązków kanałowych,
3. status pokrycia `LSM+LGR` bez legacy bridge.
