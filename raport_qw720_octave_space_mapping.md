# RAPORT QW-720: MAPOWANIE OKTAW → PRZESTRZEŃ
**Data:** 2025-12-08 23:37:28.214221

## 1. Problem
Oktawy (d=1-12) muszą być zmapowane na przestrzeń (r ~ Å) dla widma atomowego.

## 2. Testowane Hipotezy
| Hipoteza | Formuła | Średni błąd |
|----------|---------|-------------|
| Potęgowe (d^α) | - | 10194.8% |
| Eksponencjalne (exp(d/α)) | - | 9905.3% |
| Odwrotność K(d) | - | 52.8% |
| Liniowe (d×a_B) | - | 34.6% |
| **Refined (1/|K(d)|)** | **r = r₀/|K(d)|** | **52.8%** |

## 3. Finalna Funkcja Mapowania
**Liniowe (d×a_B)**

```python
def octave_to_space(d):
    # Liniowe (d×a_B)
    return mapping_func(d)
```

## 4. Weryfikacja na Promieniach Orbitali
| n | d | r_pred (Å) | r_exp (Å) | Błąd |
|---|---|------------|-----------|------|
| 1 | 1.33 | 0.5292 | 0.5292 | 0.0% |
| 2 | 9.33 | 3.7043 | 2.1167 | 75.0% |
| 3 | 17.33 | 6.8795 | 4.7626 | 44.4% |
| 4 | 25.33 | 10.0546 | 8.4668 | 18.8% |

**Średni błąd:** 34.6%

## 5. Wnioski
🟡 Mapowanie częściowo działa. Wymaga dalszych badań.
