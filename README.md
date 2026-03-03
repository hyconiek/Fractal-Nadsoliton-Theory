# Fractal Information Nadsoliton (FIN) Theory

Aktualny README projektu ToE (stan na 2026-03-03).

## 1. Czym jest ten projekt

Repozytorium zawiera rozwijaną teorię `FIN/Nadsoliton` oraz pełny pipeline badań `QW-*` (modele, audyty, bramki rygoru, raporty).

Główna idea: jedna dynamika nadsolitonu + kernel sprzężeń ma wyjaśniać emergencję obserwabli fizycznych.

## 2. Aktualny stan naukowy (uczciwy snapshot)

### Co jest mocne

- Empiryczna ścieżka confirmatory przeszła komplet bramek:
  - `QW-1852`: PASS (precheck datasetu)
  - `QW-1853`: PASS (PTA + GW)
  - `QW-1902`: `EMPIRICAL_CLOSURE_PASS` z `metric_score ~ 0.980`
- Stabilność transferowa parametru mostkującego `alpha`:
  - `QW-1912`: discovery/holdout PASS
  - `QW-1913`: multisplit PASS ALL FOLDS (`alpha=[6,6,6]`)

### Co jest nadal otwarte

- Pełne domknięcie derivacyjne ToE (bez mostu `alpha`) pozostaje niezamknięte:
  - `QW-1745`: `MICROMODEL_ITERATION_OPEN`
  - `QW-1890`: `TOE_NOT_CLOSED_REQUIRES_DERIVATIONAL_REFORMULATION`
- Wniosek: obecnie mamy silne domknięcie operacyjno-empiryczne, ale nie końcowy dowód publikacyjny ToE.

## 3. Najważniejsze artefakty (ostatnia faza)

- Plan roboczy: `PLAN_NAPRAWA_TOE_ROBOCZY.md`
- Empiryczne domknięcie:
  - `RAPORT_QW1852_EXTERNAL_CONFIRMATORY_DATA_PRECHECK.md`
  - `RAPORT_QW1853_JOINT_EXTERNAL_CONFIRMATORY_V2.md`
  - `RAPORT_QW1902_EMPIRICAL_CLOSURE_GATE.md`
- Stabilność i transfer:
  - `RAPORT_QW1912_EXTERNAL_PTA_SPLIT_VALIDATION.md`
  - `RAPORT_QW1913_EXTERNAL_PTA_MULTISPLIT_TRANSFER_STRESS.md`
- Scorecard potencjału ToE:
  - `RAPORT_QW1914_TOE_POTENTIAL_SCORECARD.md`

## 4. Szybka reprodukcja najnowszej ścieżki

Uruchamiaj w katalogu repo:

```bash
python3 QW_1910_EXTERNAL_PTA_ALPHA_ATTAINABILITY_SCAN.py
python3 QW_1911_EXTERNAL_SOURCE_DATASET_ASSEMBLY_ALPHA.py
python3 QW_1852_EXTERNAL_CONFIRMATORY_DATA_PRECHECK.py --candidate-dir external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg
python3 QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py
python3 QW_1902_EMPIRICAL_CLOSURE_GATE.py
python3 QW_1912_EXTERNAL_PTA_SPLIT_VALIDATION.py
python3 QW_1913_EXTERNAL_PTA_MULTISPLIT_TRANSFER_STRESS.py --k-folds 3
python3 QW_1914_TOE_POTENTIAL_SCORECARD.py
```

## 5. Interpretacja wyników

- Jeśli `QW-1902` = `EMPIRICAL_CLOSURE_PASS`, pipeline empiryczny jest domknięty.
- Jeśli `QW-1912/1913` utrzymują PASS, sygnał nie wygląda na jednorazowe dopasowanie splitu.
- Finalny claim ToE wymaga jeszcze niezależnej replikacji blind/preregistered na nowym pakiecie danych + domknięcia derivacyjnego bez ansatzu.

## 6. Uwaga metodologiczna

To repo zawiera bardzo dużo historycznych wersji i plików roboczych. Aktualny stan oceny teorii należy czytać przez raporty gate/audyt z końcowej fazy (`QW-1700+` i szczególnie `QW-1850+`).
