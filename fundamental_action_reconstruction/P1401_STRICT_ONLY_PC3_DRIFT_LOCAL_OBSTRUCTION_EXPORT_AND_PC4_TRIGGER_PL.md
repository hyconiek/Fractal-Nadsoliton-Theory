# P1401 — STRICT-ONLY PC3 DRIFT LOCAL OBSTRUCTION EXPORT + PC4 TRIGGER (PL)

## Cel kroku
Po `P1400` (`ROBUST_FAIL`) zamykamy klasę `PC3` uczciwie i noncyklicznie:
1. formalny eksport lokalnej przeszkody `PC3-DRIFT-v1`,
2. zamknięcie pętli `PC3` pod tym samym blocker-cut,
3. aktywacja nowej klasy provider `PC4` (strict-only).

## Zakres strict-only
- `legacy_bridge_used = false`
- brak bridge do legacy
- brak ukrytego retuningu progów po robust fail

## Wynik formalny
- `PC3_LOCAL_OBSTRUCTION_ID := PC3-DRIFT-v1`
- `PC3_LOCAL_OBSTRUCTION_STATUS := EXPORTED`
- `PC3_LOOP_STATUS := CLOSED_NONCYCLIC`
- `PC4_TRIGGER := ACTIVATED`
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Znaczenie naukowe
To domyka etap `PC3` bez fałszywego PASS i utrzymuje strategiczny rygor `QW-2381/2382/2383` (bez zapętlania tej samej rodziny poprawek).

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1402_STRICT_ONLY_PC4_PROVIDER_BASELINE_DESIGN`: zaprojektować baseline `PC4` z nową rodziną anchorów, zamrozić metryki i przygotować pierwszy run dual-metric.

## Omówienie dla laika
To jak uczciwe zakończenie wersji prototypu: zamiast w nieskończoność poprawiać wariant, który nie przechodzi testu stabilności, oficjalnie zapisujemy jego ograniczenie i przechodzimy do nowej generacji modelu.
