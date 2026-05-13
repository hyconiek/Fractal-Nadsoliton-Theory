# P1506 — S4.56 F(Nadsoliton) ⇒ L_SM + L_GR Coupled Operator Consistency Packet (PL)

Status: `P1506_EXECUTED_STRICT_COUPLED_OPERATOR_CONSISTENCY_CHECK`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok w kierunku rozwiązania `QW-2191` na torze:

`F (Nadsoliton) ⇒ L_SM + L_GR`

bez bridge do legacy, w rygorze strict-only.

## Decyzja profesorska

Po `P1505` przechodzimy z poziomu opisu na poziom testu sprzężenia:

1. sprawdzić spójność kanału `F -> LSM`,
2. sprawdzić spójność kanału `F -> LGR`,
3. sprawdzić wspólną orientację selektora (czy oba kanały mówią ten sam kierunek).

## Kryteria fizycznej spójności

Test przechodzi tylko gdy jednocześnie:

1. strict-internal selector source jest obecny i dodatni,
2. mapowanie `F_to_LSM_weight + F_to_LGR_weight = 1` (normalizacja),
3. oba kanały mają tę samą orientację selektora,
4. brak złamania guardrail: `legacy_bridge_used = false`,
5. `qw2191_closed` pozostaje `false` (bez przedwczesnego claimu).

## Wynik P1506

Publikujemy checkpoint sprzężonej spójności operatorowej jako krok
bezpośrednio przybliżający domknięcie `QW-2191` na ścieżce strict-only.
