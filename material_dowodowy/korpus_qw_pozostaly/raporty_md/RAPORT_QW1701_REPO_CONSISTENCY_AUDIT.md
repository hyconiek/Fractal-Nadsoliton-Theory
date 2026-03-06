# RAPORT QW-1701: REPO CONSISTENCY AUDIT

- Data UTC: 2026-03-02T15:09:39.388534+00:00
- Zakres: 617 skryptów `.py`, 631 dokumentów `.md`, 130 raportów `.json`
- Zakres dat plików (mtime UTC): 2025-11-12T02:03:28.361481+00:00 → 2026-03-02T15:09:14.004856+00:00

## 1) Pokrycie i mapowanie
- Skrypty z identyfikatorem QW: 458
- Skrypty zapisujące raport w kodzie: 9
- Skrypty QW z powiązanym raportem `.md/.json`: 451

## 2) Markery metodologiczne
- Skrypty z markerami optymalizacji/fittingu: 236
- Skrypty z markerami kalibracji/odniesień eksperymentalnych: 283
- Skrypty z hardcoded wartościami SM: 108

## 3) Markery raportowe
- Dokumenty typu synteza/podsumowanie/status: 459
- Dokumenty z deklaracjami `EXACT` / `0.00%`: 66
- Dokumenty z wzmianką o tautologii: 50
- Dokumenty z wzmianką o fittingu: 136

## 4) Chronologia (wg dat plików)
- Miesięczne trendy (MD):
  - 2025-11: docs=79, pos=859, neg=846, exact=62, taut=62, fit=1034
  - 2025-12: docs=520, pos=2287, neg=1630, exact=89, taut=89, fit=317
  - 2026-01: docs=18, pos=91, neg=87, exact=5, taut=0, fit=82
  - 2026-02: docs=14, pos=138, neg=132, exact=25, taut=0, fit=89
- Miesięczne trendy (PY):
  - 2025-11: py=221, opt=156, calib=166, hardcoded_sm=77
  - 2025-12: py=367, opt=75, calib=112, hardcoded_sm=31
  - 2026-02: py=28, opt=4, calib=4, hardcoded_sm=0
  - 2026-03: py=1, opt=1, calib=1, hardcoded_sm=0

## 5) Niespójności między raportami
- Liczba QW z sygnałem sprzeczności (success vs fail, lub exact vs krytyka): 725
- Przykładowe QW: QW-1, QW-2, QW-3, QW-4, QW-5, QW-6, QW-7, QW-8, QW-9, QW-10, QW-11, QW-14, QW-17, QW-18, QW-19, QW-20, QW-21, QW-22, QW-23, QW-24

## 6) Wysokie ryzyko metodologiczne (próbka)
- Skrypty z hardcoded SM (próbka):
  - `0.3 NON-TRIVIAL GROUND STATE DISCOVERY.py`
  - `0.4 IN-DEPTH ANALYSIS OF PROVIDED SOLUTION & REPORT.py`
  - `114_GENERATOR_OBSERVABLE_MAPPING.py`
  - `114_GENERATOR_OBSERVABLE_MAPPING_v2.py`
  - `115_DIAGNOSTICS.py`
  - `118_COMPOSITE_HIGGS_AND_EMERGENT_MASSES.py`
  - `122_LEPTON_MASS_UNIFIED_MECHANISM.py`
  - `16 RUNNING COUPLING CALIBRATION ANALYSIS WITH NEGATIVE BUT SCIENTIFICALLY VALUABLE RESULTS.py`
  - `19 UNIFIED GEOMETRODYNAMIC SUPERSOLITON THEORY: COMPLETE IMPLEMENTATION MAJOR BREAKTHROUGH.py`
  - `20 zadań quick win.py`
  - `3 Hierarchical Resonant Coupling for SM Mass Spectrum Reproduction.py`
  - `39mergepopr_ENHANCED.py`
  - `39mergepopr_ORIGINAL.py`
  - `41PHASE XI.2 FORMAL VERIFICATION RESULTS: CRITICAL LIMITATIONS IDENTIFIED.py`
  - `44 FOUR HIGH-PROBABILITY RESEARCH TASKS FROM FRACTAL SUPERSOLITON THEORY.py`
  - `45 IMPLEMENTATION OF UNIFIED HAMILTONIAN WITH DOUBLE-VALLEY MECHANISM.py`
  - `47 PODSUMOWANIE WYKONANIA 10 ZADAŃ KONTYNUACYJNYCH TEORII SUPERSOLITONA.py`
  - `49 MULTI-TASK ANALYSIS: ZADANIE B2, C1, D1, E1.py`
  - `5 LINKING EMERGENT GAUGE STRUCTURE TO BOSON MASSES VIA THE HIGGS MECHANISM.py`
  - `50 ZADANIE 11: Metryka z A_μ i test Poissona oraz zadanie 29.py`
- Dokumenty z jednoczesnym `EXACT` i wzmianką `tautolog/fitting` (próbka):
  - `68 ZADANIA QW-V52, QW-V53, QW-V54, QW-V55, QW-V56: ROZWIĄZANIE KRYTYCZNYCH PROBLEMÓW.md`
  - `71 ZADANIA QW-V62, QW-V63, QW-V64, QW-V65, QW-V66: FINALNE DOPRACOWANIE DO PRECYZJI.md`
  - `ANALIZA_FITTINGU_I_TRIKOW_KOMPENSACYJNYCH.md`
  - `FINAL_SUMMARY_TEORIA_WSZYSTKIEGO.md`
  - `FULL_LOG_COMPRESSED_EXTREME_QW420_1200.md`
  - `FULL_LOG_COMPRESSED_PHASE1_AUDIT.md`
  - `FULL_LOG_COMPRESSED_PHASE2_UPGRADE.md`
  - `FULL_LOG_OPERATIONAL_PHASE2.md`
  - `FULL_LOG_SPINORS_PHASE1.md`
  - `INDEKS_PROJEKTU_BADAN_1_118_FINALNA.md`
  - `KONTEXT_TEORII_DLA_AI_RESEARCH.md`
  - `KRYTYKA_BRAKI_NIEZGODNOSCI.md`
  - `OSIAGNIECIE_HISTORYCZNE_TEORIA_WSZYSTKO.md`
  - `OSTATECZNA_WIADOMOSC.md`
  - `OSTATECZNY_RAPORT_14_XI_2025.md`
  - `PODSUMOWANIE_BADAN_116_117_118.md`
  - `QW-335-339_FINAL_VERDICT.md`
  - `RAPORT_AUDYT_METODOLOGICZNY_QW500_QW826.md`
  - `RAPORT_RYGORYSTYCZNEJ_WERYFIKACJI.md`
  - `TEMP_FULL_REPORT_CONTEXT.md`

## 7) Artefakty
- JSON szczegółowy: `report_qw1701_repo_consistency_audit.json`

## Wniosek (techniczny)
Repo zawiera równolegle silne deklaracje potwierdzeń i liczne raporty krytyczne/falsyfikacyjne. W praktyce status teorii jest heterogeniczny i wymaga wersjonowanego rejestru hipotez z rozróżnieniem: wynik predykcyjny vs wynik dopasowany/diagnostyczny.
