# QW-1613: Mass Generations from T(21,3) Deformation

**Data:** 2025-12-29 22:27:44

**Status:** FAILED

## Model
m(d) = m_e × exp(-γ × (d - 6.0))

gdzie γ = 1.52

## Wyniki

| Stosunek | Predykcja | Eksperyment | Błąd |
|----------|-----------|-------------|------|
| m_μ/m_e | 44.70 | 206.77 | 78.4% |
| m_τ/m_e | 298.87 | 3477.22 | 91.4% |
| m_τ/m_μ | 6.69 | 16.82 | 60.2% |

## Raw Log
```
================================================================================
QW-1613: MASS GENERATIONS FROM T(21,3) DEFORMATION
================================================================================
Data: 2025-12-29 22:27:44

[1] WĘZŁY TORUSOWE DLA LEPTONÓW
------------------------------------------------------------
Elektron: T(21,3), Q = 24, d = 6.0
Mion:     T(8,5), Q = 13, d = 3.5
Tau:      T(5,3), Q = 8, d = 2.25

[2] PRZEWIDYWANE MASY (KALIBRACJA: m_e)
------------------------------------------------------------
m_e (kalibracja):  0.5110 MeV
m_μ (predykcja):   22.8423 MeV
m_τ (predykcja):   152.7212 MeV

[3] STOSUNKI MAS
------------------------------------------------------------
Stosunek           Predykcja  Eksperyment     Błąd %
--------------------------------------------------
m_μ/m_e                44.70       206.77      78.38%
m_τ/m_e               298.87      3477.22      91.40%
m_τ/m_μ                 6.69        16.82      60.24%

[4] WERDYKT KOŃCOWY
============================================================
❌ Średni błąd stosunków: 76.7%
   Model wymaga fundamentalnej rewizji

UWAGA: To jest 'phenomenological consistency check', nie dowód TOE.
OVERALL STATUS: FAILED
```
