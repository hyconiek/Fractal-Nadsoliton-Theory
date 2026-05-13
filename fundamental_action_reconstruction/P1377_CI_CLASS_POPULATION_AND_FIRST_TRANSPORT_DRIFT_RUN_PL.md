# P1377 c_i Class Population And First Transport-Drift Run (PL)

Status: `P1377_EXECUTED_FIRST_DRIFT_RUN_INCOMPLETE_CLASS_NO_FALSE_PASS`
As of: `2026-05-12`

## Cel

Po `P1376` wykonujemy pierwszy uruchamialny run metryki `transport_drift`
z częściową populacją klasy `c_i`.

To jest krok pośredni przed theorem-level werdyktem `L-B1-02`.

## Dane robocze (v1 trial)

Użyte współczynniki:

- `c_g3`: A=1.00, B=1.03
- `c_g2`: A=1.00, B=0.98
- `c_g1`: A=1.00, B=1.01

Brak jeszcze `c_mix` -> klasa niepełna.

## Wynik metryki

`transport_drift = max(0.03, 0.02, 0.01) = 0.03`

`epsilon_transport` nadal `UNSET` (zgodnie z `P1375/P1376`),
więc nie wolno wydać werdyktu PASS/FAIL theorem-level.

## Werdykt

`L-B1-02_STATUS := INCOMPLETE_CLASS_TRIAL_ONLY`

Znaczenie:

1. pipeline drift działa i jest liczony jawnie,
2. ale brak `c_mix` i brak zatwierdzonego progu `epsilon_transport`
   blokuje końcowy wniosek o inwariancji.

## Decyzja profesorska

Następny krok MUSI domknąć:

1. pełną populację klasy (`+ c_mix`),
2. kalibrację i freeze `epsilon_transport`,
3. dopiero potem pierwszy formalny PASS/FAIL `L-B1-02`.

