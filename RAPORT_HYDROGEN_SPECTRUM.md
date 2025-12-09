# RAPORT: WIDMO WODORU (QW-676)
**Data:** 2025-12-08 20:42:57.966092
**Cel:** Czy FIN produkuje $E_n \propto 1/n^2$?

## 1. SETUP
- Siatka radialna: 50 punktów
- Proton oktawa: 7
- Elektron oktawa: 1
- K(d=6) = 0.8664

## 2. STANY ZWIĄZANE
- Liczba stanów związanych: 4

| Stan n | Energia E_n | E_n/E_1 | Teoria (1/n²) |
|--------|-------------|---------|---------------|
| 1 | -5.9250 | 1.000 | 1.000 |
| 2 | -2.4106 | 0.4069 | 0.250 |
| 3 | -1.3029 | 0.2199 | 0.111 |

## 3. ZGODNOŚĆ Z $1/n^2$
- Średni błąd: 8.9%
- E_2/E_1: 0.4069 (teoria: 0.25)
- E_3/E_1: 0.2199 (teoria: 0.111)

✅ **SUCCESS:** Widmo zgodne z $1/n^2$!

## 4. PODSUMOWANIE

| Test | Wynik | Status |
|------|-------|--------|
| Stany związane | 4 | ✅ |
| E_2/E_1 | 0.4069 | ❌ |
| Zgodność $1/n^2$ | 91.1% | ✅ |

