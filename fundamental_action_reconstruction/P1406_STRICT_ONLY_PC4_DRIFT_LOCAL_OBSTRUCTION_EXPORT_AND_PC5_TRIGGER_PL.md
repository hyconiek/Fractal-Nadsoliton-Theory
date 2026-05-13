# P1406 — STRICT-ONLY PC4 DRIFT LOCAL OBSTRUCTION EXPORT + PC5 TRIGGER (PL)

## Cel kroku
Po `P1405` (`ROBUST_FAIL`) zamykamy pętlę `PC4` w trybie noncyklicznym:
1. formalny eksport przeszkody lokalnej `PC4-DRIFT-v1`,
2. zamknięcie klasy poprawek `PC4`,
3. aktywacja nowej klasy provider `PC5` bez bridge do legacy.

## Zakres strict-only
- `legacy_bridge_used = false`
- `noncyclic_provider_shift = enforced`
- brak retuningu progów po robust fail

## Wynik formalny
- `PC4_LOCAL_OBSTRUCTION_ID := PC4-DRIFT-v1`
- `PC4_LOCAL_OBSTRUCTION_STATUS := EXPORTED`
- `PC4_LOOP_STATUS := CLOSED_NONCYCLIC`
- `PC5_TRIGGER := ACTIVATED`
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1407_STRICT_ONLY_PC5_PROVIDER_BASELINE_DESIGN` i zamrozić nową bazę `PC5` przed pierwszym runem dual-metric.

## Omówienie dla laika
To uczciwe zakończenie obecnej wersji testu: skoro wynik nie jest stabilny, zapisujemy ograniczenie i przechodzimy do nowej generacji modelu, zamiast poprawiać bez końca ten sam wariant.
