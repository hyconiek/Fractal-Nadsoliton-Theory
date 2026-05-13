# P1511 — S4.61 Admissible-Perturbation Boundary Map F⇒LSM+LGR Packet (PL)

Status: `P1511_EXECUTED_ADMISSIBLE_BOUNDARY_MAP_STRICT_ONLY`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok po `P1510`: wyznaczyć mapę granicy perturbacji
na torze `F(Nadsoliton)=>LSM+LGR`, która rozdziela:

1. strefę robust (draft theorem utrzymany),
2. strefę odrzucenia (gałąź kontradykcji aktywna).

## Decyzja profesorska

Przechodzimy z jakościowego opisu do ilościowego mapowania granicy. To jest
warunek rygoru fizycznego przed jakimkolwiek ruchem statusowym.

## Zakres strict-only

Mapa używa wyłącznie obiektów strict-side (`S_internal_v1`, `W_Fmap_v1`) i
nie przenosi żadnego legacy bridge claimu.

## Kryterium

Dla każdej próby perturbacji oceniamy A1..A5 z `P1508`.

- `ALL_TRUE` => punkt robust,
- `ANY_FALSE` => punkt odrzucenia.

## Wynik P1511

Publikujemy liczbową granicę robust/rejection jako podstawę do następnego kroku:
koncentracji na minimalnym zestawie warunków potrzebnym do trwałego
`F(Nadsoliton)=>LSM+LGR`.
