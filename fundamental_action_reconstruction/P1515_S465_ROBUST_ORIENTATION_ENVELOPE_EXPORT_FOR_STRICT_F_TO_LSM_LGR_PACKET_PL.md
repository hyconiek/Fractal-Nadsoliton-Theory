# P1515 — S4.65 Robust-Orientation Envelope Export For Strict F⇒LSM+LGR Packet (PL)

Status: `P1515_EXECUTED_ROBUST_ORIENTATION_ENVELOPE_EXPORT`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok po `P1514`: wyeksportować formalną obwiednię
orientacyjną (robust-orientation envelope) dla strict-side draftu
`F(Nadsoliton)=>LSM+LGR`.

## Decyzja profesorska

Stabilne orientacje muszą być jawnie wylistowane jako aneks dopuszczalności,
aby każdy dalszy claim był kontrolowalny i falsyfikowalny.

## Definicja robocza envelope

Do obwiedni wchodzą tylko stany, które przechodzą `A4*`:

`shared_orientation AND orientation_margin >= 0.5`.

## Wynik P1515

Publikujemy obwiednię orientacji robust oraz listę stanów poza obwiednią
jako formalny załącznik do strict coupled theorem draft.
