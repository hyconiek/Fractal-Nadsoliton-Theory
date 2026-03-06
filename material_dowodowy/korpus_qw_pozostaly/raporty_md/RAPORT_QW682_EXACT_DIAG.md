# RAPORT: QW-682 EXACT DIAGONALIZATION BELL TEST
**Data:** 2025-12-06 23:17:49.323832
**Cel:** Test Bell inequality z PEŁNĄ mechaniką kwantową (nie mean-field)

## 1. METODOLOGIA

- **N spins:** 8 (Hilbert space dim = 256)
- **Topology:** Triangular (frustrated)
- **Coupling:** J = -1.0 (antiferromagnetic)
- **Method:** Full exact diagonalization (scipy.linalg.eigh)

## 2. GROUND STATE

- **Energy:** E_0 = -3.000000
- **Gap:** ΔE = 0.000000
- **Entanglement entropy:** S_vN = 1.1033

## 3. CHSH TEST

| Para | S (CHSH) | Violation? |
|------|----------|------------|
| (0, 4) | 0.6561 | ❌ |
| Best (1, 7) | 0.6561 | ❌ |
| Mean | 0.6561 | - |

**Classical bound:** S ≤ 2.0
**Tsirelson bound:** S ≤ 2√2 ≈ 2.828

## 4. KORELACJE (para 0-4)

- E(a,b) = -0.050656
- E(a,b') = 0.050360
- E(a',b) = 0.378405
- E(a',b') = 0.378701

## 5. PODSUMOWANIE

### ❌ PORAŻKA: NO BELL VIOLATION

Max S = 0.6561 ≤ 2.0

**Przyczyna:** Ground state nie jest maksymalnie splątany.
Frustracja tworzy spin liquid, ale nie singlet-like correlations.

**Rekomendacja:**
1. Zwiększyć N (więcej sprzężeń)
2. Próbować innych topologii (kagome, triangular 2D)
3. Dodać anisotropię (XXZ model)

## 6. PORÓWNANIE Z POPRZEDNIMI

| Badanie | Metoda | Max S | Status |
|---------|--------|-------|--------|
| QW-680 | Mean-field | 0.17 | ❌ |
| **QW-682** | **Exact Diag** | **0.66** | ❌ |
