# QW-1617/1618/1619: Optional Studies (Combined)

**Data:** 2025-12-29 22:27:44

## Summary

| Study | Topic | Status |
|-------|-------|--------|
| QW-1617 | Running α(Q²) | PARTIAL |
| QW-1618 | Preon Stability | VERIFIED |
| QW-1619 | CKM from Knots | PARTIAL |

## Notes
- These studies are **model-dependent** and **speculative**
- They serve as **consistency checks**, not proofs
- QW-1619 carries **high numerology risk**

## Raw Log
```
================================================================================
QW-1617/1618/1619: OPTIONAL STUDIES (COMBINED)
================================================================================
Data: 2025-12-29 22:27:44

============================================================
QW-1617: RUNNING COUPLING α(Q²)
============================================================
   Q [MeV]       α(QED)       α(FIN)     Diff %
--------------------------------------------------
       0.5     0.007297     0.007286       0.15%
       1.0     0.007305     0.007292       0.18%
      10.0     0.007331     0.007312       0.27%
      91.2     0.007356     0.007330       0.36%
    1000.0     0.007384     0.007351       0.45%

QW-1617 STATUS: PARTIAL (qualitative consistency only)

============================================================
QW-1618: PREON STABILITY
============================================================
  n_preons      E_total          E/n    Stable?
--------------------------------------------------
         1       1.1000       1.1000          ✅
         2      -0.3726      -0.1863          ✅
         3      -4.4178      -1.4726          ✅
         4     -11.0355      -2.7589          ✅
         5     -20.2259      -4.0452          ✅
         6     -31.9888      -5.3315          ✅
         7     -46.3244      -6.6178          ✅

QW-1618 STATUS: VERIFIED

============================================================
QW-1619: CKM FROM KNOT ROTATIONS
============================================================
Przybliżona struktura CKM (qualitative):

|V_ud| ≈ 0.8187 (exp: 0.974)
|V_us| ≈ 0.2019 (exp: 0.225)
|V_cd| ≈ 0.2466 (exp: 0.221)
|V_cs| ≈ 1.0000 (exp: 0.973)

Hierarchia |V_ud| > |V_us|: ✅
QW-1619 STATUS: PARTIAL (qualitative mapping only)
UWAGA: Wysokie ryzyko numerologii!

============================================================
PODSUMOWANIE OPTIONAL STUDIES
============================================================
QW-1617 (Running α):    PARTIAL
QW-1618 (Preon Stab.):  VERIFIED
QW-1619 (CKM Knots):    PARTIAL
```
