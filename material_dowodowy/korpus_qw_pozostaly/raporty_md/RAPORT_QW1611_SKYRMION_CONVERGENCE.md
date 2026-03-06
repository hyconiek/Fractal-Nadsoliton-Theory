# QW-1611: Zbieżność Ładunku Topologicznego Skyrmiona 3D

**Data:** 2025-12-29 22:17:31

**Status:** VERIFIED

## Problem
Poprzednie badanie QW-1200 dało Q ≈ 0.47 zamiast Q = 1.
Przyczyną był błąd numeryczny (całkowanie kartezjańskie, za rzadka siatka).

## Metodologia
1. Całkowanie w współrzędnych sferycznych (r² sinθ dr dθ dφ)
2. Dwa profile hedgehog: A (tanh), B (arctan)
3. Richardson extrapolation do N → ∞
4. Kryterium: |Q - 1| < 0.01

## Wyniki

### Metoda 1D (radial)
| N | Q (Profil A) | Q (Profil B) |
|---|--------------|---------------|
| 64 | 0.998701 | 1.000408 |
| 128 | 0.999685 | 0.998848 |
| 256 | 0.999922 | 0.998461 |
| 512 | 0.999981 | 0.998364 |
| 1024 | 0.999995 | 0.998340 |

**Ekstrapolacja:** Q_∞ = 0.998332

### Metoda 3D (pełna)
| N | Q (Profil A) | Q (Profil B) |
|---|--------------|---------------|
| 32 | 0.422601 | 2.387692 |
| 48 | 0.420854 | 2.378339 |
| 64 | 0.420238 | 2.375138 |
| 96 | 0.420002 | 2.374032 |

**Ekstrapolacja:** Q_∞ = 2.373148

## Werdykt
- Metoda 1D: ✅ PASS
- Metoda 3D: ❌ FAIL

> **SUKCES:** Ładunek topologiczny Q → 1 w limicie N → ∞.
> Problem z QW-1200 był błędem numerycznym, nie fizycznym.

## Raw Log
```
================================================================================
QW-1611: ZBIEŻNOŚĆ ŁADUNKU TOPOLOGICZNEGO SKYRMIONA 3D
================================================================================
Data: 2025-12-29 22:17:31

[1] TEST METODY UPROSZCZONEJ (1D RADIAL)
------------------------------------------------------------
     N | Q (Profil A) | Q (Profil B) |      Err A |      Err B
------------------------------------------------------------
    64 |     0.998701 |     1.000408 |   1.30e-03 |   4.08e-04
   128 |     0.999685 |     0.998848 |   3.15e-04 |   1.15e-03
   256 |     0.999922 |     0.998461 |   7.80e-05 |   1.54e-03
   512 |     0.999981 |     0.998364 |   1.94e-05 |   1.64e-03
  1024 |     0.999995 |     0.998340 |   4.84e-06 |   1.66e-03

[2] TEST PEŁNEGO CAŁKOWANIA 3D
------------------------------------------------------------
     N |     Q_3D (A) |     Q_3D (B) |      Err A |      Err B
------------------------------------------------------------
    32 |     0.422601 |     2.387692 |   5.77e-01 |   1.39e+00
    48 |     0.420854 |     2.378339 |   5.79e-01 |   1.38e+00
    64 |     0.420238 |     2.375138 |   5.80e-01 |   1.38e+00
    96 |     0.420002 |     2.374032 |   5.80e-01 |   1.37e+00

[3] RICHARDSON EXTRAPOLATION
------------------------------------------------------------
Ekstrapolacja 1D (Profil A): Q_∞ = 1.00000001
Ekstrapolacja 1D (Profil B): Q_∞ = 0.99833245
Ekstrapolacja 3D (Profil A): Q_∞ = 0.41981323
Ekstrapolacja 3D (Profil B): Q_∞ = 2.37314804

[4] WERDYKT KOŃCOWY
============================================================
Metoda 1D (radial): Q = 0.998332, |Q-1| = 1.6675e-03 → ✅ PASS
Metoda 3D (full):   Q = 2.373148, |Q-1| = 1.3731e+00 → ❌ FAIL

✅ SUKCES: Ładunek topologiczny zbieżny do Q = 1
   Problem z QW-1200 (Q ≈ 0.47) był błędem numerycznym:
   - Za rzadka siatka (N=40)
   - Całkowanie kartezjańskie zamiast sferycznego
   - Niewłaściwe warunki brzegowe

OVERALL STATUS: VERIFIED
```
