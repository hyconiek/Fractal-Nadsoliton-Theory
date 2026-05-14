# P1665 / S615 Spin-2 Propagator Unitarity Proxy Packet (Strict-only)

Status: `P1665_EXECUTED_SPIN2_PROPAGATOR_PROXY_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel

Wykonać kolejny uczciwy krok po `P1664`: z pełnego strictowego `L_GR` wyprowadzić
proxy-analizę propagatora spin-2 i jawnie zaznaczyć, co jest jeszcze brakującym theoremem.

## Wejście z toru strict

`K_strict -> (Mpl2, cR2, cRic2, cRiem2) -> L_GR`

gdzie:

`L_GR = sqrt(-g)[(Mpl2/2)R + cR2 R^2 + cRic2 Ricci^2 + cRiem2 Riemann^2]`

## Analiza proxy spin-2

Wokół tła Minkowskiego używamy standardowego rozszczepienia projektorowego (operacyjny poziom):

- biegun bezmasowy spin-2: `k^2 = 0`,
- ciężki biegun spin-2: `m2_2 ~ Mpl2 / (cRic2 + 4*cRiem2)`.

Proxy-warunki:

1. `Mpl2 > 0` (prawidłowa dodatnia energia gałęzi bezmasowej),
2. `m2_2 > 0` (brak tachionu w gałęzi ciężkiej),
3. znak residuum ciężkiej gałęzi zaznaczony jawnie jako `ghost_risk_flag`.

## Zasada bez fałszywego pass

Nawet jeśli 1-2 przejdą, brak pełnego theoremu unitarności operatorowej
=> status globalny pozostaje `OPEN_OBLIGATION`.

## Rekomendowany następny uczciwy krok

`S616`: policzyć pełny operator propagatora dla klasy teł FRW/stałej krzywizny i dodać test spektralny (bieguny + residua) dla całego sektora spin-2/spin-0.

## Omówienie dla laika

To jak sprawdzanie, czy w układzie drgań nie pojawia się „fałszywy rezonans”,
który dawałby niefizyczną energię. Tu robimy pierwszy techniczny test tej stabilności,
ale pełny certyfikat bezpieczeństwa nadal wymaga dalszych dowodów.
