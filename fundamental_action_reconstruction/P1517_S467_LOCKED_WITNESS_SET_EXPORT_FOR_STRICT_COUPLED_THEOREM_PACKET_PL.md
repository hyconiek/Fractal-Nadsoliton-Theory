# P1517 — S4.67 Locked Witness Set Export For Strict Coupled Theorem Packet (PL)

Status: `P1517_EXECUTED_LOCKED_WITNESS_SET_EXPORT`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok po `P1516`: uruchomić rerun sprzężonego
strict-theorem z aktywnym lockiem aneksu i wyeksportować `locked witness set`.

## Decyzja profesorska

Wszystkie oceny twierdzenia muszą przejść przez lock (`robust_envelope`).
Niezależnie od wyniku lokalnego, brak claimu globalnego i brak bridge do legacy.

## Wynik P1517

Publikujemy zbiór witnessów dopuszczonych przez lock i zbiór odrzuconych
scenariuszy jako formalny wynik rerunu `F(Nadsoliton)=>LSM+LGR`.
