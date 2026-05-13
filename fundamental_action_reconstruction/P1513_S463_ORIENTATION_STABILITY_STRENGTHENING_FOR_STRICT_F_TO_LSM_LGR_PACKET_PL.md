# P1513 — S4.63 Orientation-Stability Strengthening For Strict F⇒LSM+LGR Packet (PL)

Status: `P1513_EXECUTED_ORIENTATION_STABILITY_STRENGTHENING`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok po `P1512`: wzmocnić krytyczny warunek
`A4_shared_orientation` na torze `F(Nadsoliton)=>LSM+LGR`.

## Decyzja profesorska

Skoro `A4` jest aktywnym warunkiem granicznym odrzucenia, wzmacniamy dowód
przez jawny margines stabilności orientacji zamiast rozszerzać claim.

## Metoda

1. Definiujemy orientacyjne perturbacje kanałów (`SM_preferred`, `GR_preferred`,
   `mixed`, `neutral`),
2. testujemy utrzymanie A1..A5,
3. wyznaczamy minimalny margines orientacyjny, w którym twierdzenie pozostaje
   dopuszczone.

## Wynik P1513

Publikujemy tabelę stabilności orientacji i formalny warunek wzmocniony:

`A4* : shared_orientation AND orientation_margin >= m_min`.
