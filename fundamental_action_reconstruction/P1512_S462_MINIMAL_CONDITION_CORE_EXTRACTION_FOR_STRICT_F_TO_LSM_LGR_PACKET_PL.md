# P1512 — S4.62 Minimal-Condition Core Extraction For Strict F⇒LSM+LGR Packet (PL)

Status: `P1512_EXECUTED_MINIMAL_CORE_EXTRACTION_STRICT_ONLY`
As of: `2026-05-13`

## Cel

Wykonać kolejny uczciwy krok po `P1511`: wyznaczyć minimalny rdzeń warunków,
który podtrzymuje robust strefę dla `F(Nadsoliton)=>LSM+LGR`.

## Decyzja profesorska

Zamiast rozszerzać model, redukujemy warunki do najmniejszego zestawu
koniecznego do utrzymania spójności strict-side.

## Metoda

Używamy mapy z `P1511` i porównujemy punkty robust/rejection, by określić:

1. warunki konieczne,
2. warunki techniczne (pomocnicze),
3. jedyny krytyczny warunek brzegowy aktywujący odrzucenie.

## Wynik P1512

Publikujemy minimal core A* dla aktualnego theorem draftu oraz wskazujemy
następny krok: wzmocnienie dowodu wokół krytycznego warunku orientacji.
