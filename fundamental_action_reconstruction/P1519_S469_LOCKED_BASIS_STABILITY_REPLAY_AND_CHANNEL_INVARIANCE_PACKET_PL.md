# P1519 — S4.69 Locked-Basis Stability Replay And Channel-Invariance Packet (PL)

Status: `P1519_EXECUTED_LOCKED_BASIS_STABILITY_REPLAY`
As of: `2026-05-13`

## Cel

Wykonać kolejny uczciwy krok po `P1518`: przeprowadzić replay stabilności
zablokowanej bazy witnessów na rozszerzonej siatce perturbacji i sprawdzić
niezmienność pokrycia kanałów `F->LSM` i `F->LGR`.

## Decyzja profesorska

Utrzymujemy pełny rygor strict-only: tylko przypadki w lock-bazie,
bez legacy bridge i bez podnoszenia claimu global closure.

## Wynik P1519

Publikujemy raport stabilności bazy witnessów oraz werdykt, czy pokrycie
kanałów jest niezmienne pod testowanymi perturbacjami.
