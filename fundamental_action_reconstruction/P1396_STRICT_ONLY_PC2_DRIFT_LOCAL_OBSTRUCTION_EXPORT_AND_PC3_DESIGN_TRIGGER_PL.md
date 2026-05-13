# P1396 — STRICT-ONLY PC2 DRIFT LOCAL OBSTRUCTION EXPORT + PC3 DESIGN TRIGGER (PL)

## Cel kroku
Po `P1395` (robust fail driftu selektora w klasie `PC2`) wykonujemy **uczciwe zamknięcie lokalnej klasy poprawek**:
1. formalny eksport przeszkody lokalnej `PC2-DRIFT-v1`,
2. uruchomienie nowej, niecyklicznej gałęzi `PC3` (bez bridge do legacy),
3. utrzymanie rygoru `no-false-pass` i guardrail `QW-2381/2382/2383`.

## Zakres strict-only
- `legacy_bridge_used = false`
- brak transferu ról legacy -> strict
- brak nowych aksjomatów maskujących blocker drift

## Dane wejściowe
- `P1395_STRICT_ONLY_PC2_DRIFT_EPSILON_EDGE_ROBUSTNESS_RUN_PL.md`
- `generated/p1395_strict_only_pc2_drift_epsilon_edge_robustness_run_summary.json`

## Wynik formalny
- `PC2_LOCAL_OBSTRUCTION_ID := PC2-DRIFT-v1`
- `PC2_LOCAL_OBSTRUCTION_STATUS := EXPORTED`
- `PC2_LOOP_STATUS := CLOSED_NONCYCLIC`
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Znaczenie naukowe
To nie jest krok „wstecz”, tylko krok jakościowy: zamiast kręcić tę samą pętlę PC2 przy aktywnym robust fail, jawnie eksportujemy obstruction i otwieramy nową klasę provider/anchor (`PC3`) z czystym baseline. To realizuje noncyclic discipline i zwiększa szansę realnego postępu theorem-level.

## Trigger dla PC3
- `PC3_TRIGGER := ACTIVATED`
- `PC3_POLICY := NEW_PROVIDER_CLASS_REQUIRED`
- `PC3_FIRST_PACKET := P1397_STRICT_ONLY_PC3_PROVIDER_BASELINE_DESIGN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1397_STRICT_ONLY_PC3_PROVIDER_BASELINE_DESIGN`: zdefiniować nową klasę providerów/anchorów niezależną od rodziny PC2, zamrozić metryki (`sign_flip_rate`, `selector_drift`) i dopiero po tym wykonać pierwszy run na PC3.

## Omówienie dla laika
To jak uczciwa decyzja inżynierska: skoro jedna konstrukcja regularnie nie przechodzi testu odporności, nie „dokręcamy śrubki bez końca”, tylko oficjalnie zapisujemy ograniczenie i projektujemy nową wersję od podstaw.
