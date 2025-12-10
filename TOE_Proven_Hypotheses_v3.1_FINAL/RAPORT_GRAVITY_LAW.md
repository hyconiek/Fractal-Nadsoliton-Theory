# RAPORT: PRAWO GRAWITACJI (QW-675)
**Data:** 2025-12-06 22:09:16.684066
**Cel:** Czy FIN produkuje $F \propto 1/r^2$ (Newton)?

## 1. SETUP
- Węzły: 400
- Oktawy: 0-11
- $\alpha_{geo} = 2.7726$

## 2. MASA CENTRALNA
- Pozycja: [0, 0, 0]
- Promień: 2.0
- Węzły masowe: 7

## 3. POMIARY
| Odległość r | Energia E | Siła F=-dE/dr |
|-------------|-----------|---------------|
| 3.0 | 12.94 | -0.0000 |
| 4.0 | 12.94 | 9.1592 |
| 5.0 | -5.38 | 4.2233 |
| 6.0 | 4.49 | -3.3268 |
| 7.0 | 1.27 | 1.6091 |
| 8.0 | 1.27 | 0.2174 |
| 10.0 | -0.03 | 0.3261 |
| 12.0 | -0.03 | -0.0000 |

## 4. DOPASOWANIE PRAWA POTĘGOWEGO
**Wynik:** $F = 1275.2404 \cdot r^{-3.56}$

- Wykładnik n = **-3.56**
- Oczekiwany (Newton): n = **-2.0**
- Błąd: 1.56

🟡 **PARTIAL:** Przyciąganie, ale n=-3.56 ≠ -2.0

## 5. ALTERNATYWNA ANALIZA: E(r)
**Wynik:** $E = 62.6724 \cdot r^{-1.08}$

- Wykładnik m = **-1.08**
- Oczekiwany (Newton): m = **-1.0**

## 6. PODSUMOWANIE

| Test | Wynik | Status |
|------|-------|--------|
| Wykładnik n (siła) | -3.56 | ❌ |
| Kierunek siły | Przyciąganie | ✅ |

**Wniosek:**
FIN Theory NIE produkuje standardowego prawa grawitacji.
Możliwe przyczyny:
1. Sieć dyskretna nie ma ciągłego limitu
2. Jądro K(d) dominuje nad odległością przestrzenną
3. Potrzebna większa sieć (N > 10⁴)
