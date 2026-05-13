# P1378 c_mix Population And epsilon_transport Freeze (PL)

Status: `P1378_EXECUTED_FREEZE_PREPASS_NO_FALSE_PASS`
As of: `2026-05-12`

## Cel

Domknąć brakujący etap po `P1377`:

1. dopełnić klasę `c_i` o `c_mix`,
2. zamrozić roboczy próg `epsilon_transport` dla pierwszego pełnego testu,
3. przygotować formalny run `PASS/FAIL` w następnym pakiecie.

## Populacja c_mix (trial v1)

Ustalono roboczo:

- `c_mix(A) = 1.00`
- `c_mix(B) = 0.97`

co daje `drift_mix = 0.03`.

## Freeze progu

Dla pierwszego pełnego runu przyjmujemy:

`epsilon_transport_v1 = 0.035`

Uwaga rygorystyczna: to jest freeze operacyjny v1, a nie finalna stała
fizyczna programu.

## Werdykt tego pakietu

`L-B1-02_READINESS := FULL_CLASS_AND_EPSILON_FROZEN`

To nadal NIE jest końcowy PASS/FAIL inwariancji — tylko gotowość do niego.

## Decyzja profesorska

Następny krok: wykonać pełny formalny run `L-B1-02` na kompletnej klasie
z `epsilon_transport_v1`, bez żadnych dodatkowych korekt ad hoc.

