# Raport 112 — Deep Nadsoliton Probe
**Autor:** Krzysztof Żuchowski

Data: 2025-11-14T20:43:37.785900+00:00

## Task 1 — Algebraic probe (N=24)
Wnioski po chłopsku:
 - Słaba zgodność z prostym zamknięciem na bazie projekcji top-4 → prawdopodobnie potrzebna większa baza (więcej trybów) lub inna baza generatorów.
Metryka: fraction_below_1e-2 = 0.250, fraction_below_1e-1 = 0.250

## Task 2 — Participation Ratio scaling
Wnioski po chłopsku:
 - mode_1: PR ~ N^0.998  (jeśli alpha≈1 → rozciągły tryb; alpha≈0 → lokalny)
 - mode_2: PR ~ N^0.998  (jeśli alpha≈1 → rozciągły tryb; alpha≈0 → lokalny)
 - mode_3: PR ~ N^0.983  (jeśli alpha≈1 → rozciągły tryb; alpha≈0 → lokalny)
 - mode_4: PR ~ N^0.976  (jeśli alpha≈1 → rozciągły tryb; alpha≈0 → lokalny)
Szczegóły: patrz sekcja danych w JSON

## Task 3 — Topological defect probe
 - N=12: PR_orig=8.007, PR_def=6.235, delta_lambda=2.731903
 - N=16: PR_orig=10.671, PR_def=8.896, delta_lambda=2.687378
 - N=20: PR_orig=13.334, PR_def=11.543, delta_lambda=2.696466
 - N=24: PR_orig=15.994, PR_def=14.204, delta_lambda=2.654721
 - N=28: PR_orig=18.653, PR_def=16.851, delta_lambda=2.660711
 - N=32: PR_orig=21.310, PR_def=19.509, delta_lambda=2.621453
Wnioski: jeśli PR_def << PR_orig i delta_lambda znaczące → lokalizowane jądro topologiczne.

## Task 4 — Generator reconstruction (N=24)
 - Effective rank (liczba niezależnych generatrów w przestrzeni S(phi)) = 2
Wnioski: mała liczba generatorów sugeruje niskowymiarową algebrę operatorów; duża liczba → bardziej złożona struktura.

## Task 5 — RG landscape sweep
 - N=12: sign_changes in dg/dlns = 0
 - N=16: sign_changes in dg/dlns = 0
 - N=20: sign_changes in dg/dlns = 0
 - N=24: sign_changes in dg/dlns = 0
 - N=28: sign_changes in dg/dlns = 0
 - N=32: sign_changes in dg/dlns = 0
Wnioski: regiony z change sign wskazują na przejścia screening/antiscreening; jeśli istnieją dla pewnych N→faza.

---
Pliki wygenerowane: JSON oraz MD.
