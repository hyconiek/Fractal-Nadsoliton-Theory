# P1410 — STRICT-ONLY PC5 DRIFT LOCAL OBSTRUCTION EXPORT + PC6 TRIGGER (PL)

## Cel kroku
Po `P1409` (`ROBUST_FAIL`) formalnie zamykamy pętlę `PC5`:
1. eksport przeszkody lokalnej `PC5-DRIFT-v1`,
2. zamknięcie `PC5` noncyklicznie,
3. aktywacja nowej klasy provider `PC6` (strict-only).

## Zakres strict-only
- `legacy_bridge_used = false`
- brak bridge do legacy
- brak retuningu progów po robust-fail

## Wynik formalny
- `PC5_LOCAL_OBSTRUCTION_ID := PC5-DRIFT-v1`
- `PC5_LOCAL_OBSTRUCTION_STATUS := EXPORTED`
- `PC5_LOOP_STATUS := CLOSED_NONCYCLIC`
- `PC6_TRIGGER := ACTIVATED`
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1411_STRICT_ONLY_PC6_PROVIDER_BASELINE_DESIGN` i zamrozić nową bazę `PC6` przed pierwszym runem dual-metric.

## Omówienie dla laika
To uczciwe zamknięcie etapu: skoro obecna wersja nie przeszła testu odporności, zapisujemy to formalnie i przechodzimy do nowej wersji modelu zamiast kręcić się w tej samej pętli.
