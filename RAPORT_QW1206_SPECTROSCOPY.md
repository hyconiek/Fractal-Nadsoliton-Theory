# QW-1206: Spektroskopia Węzłów Torusowych

**Data:** 2025-12-11 02:10:03.521505

```
================================================================================
QW-1206: SPEKTROSKOPIA TOPOLOGICZNA WĘZŁÓW TORUSOWYCH
================================================================================

[1] BADANIE WIDMA DLA RÓŻNYCH WĘZŁÓW (DLA Q=24)
--------------------------------------------------------------------------------
Węzeł      Q     Score (MSE)     Slope      Widmo (pierwsze 4 mode)                 
--------------------------------------------------------------------------------
T(21,3) - Pominięto (nie jest węzłem, gcd=3)
T(13,11)  24    0.005966        0.13       [1.   1.   1.22 1.22]                    
T(19,5)   24    0.002056        0.10       [1.   1.16 1.16 1.36]                    
T(17,7)   24    0.001914        0.10       [1.   1.   1.19 1.19]                    
T(12,12) - Pominięto (nie jest węzłem, gcd=12)
T(13,8)   21    0.003680        0.14       [1.   1.25 1.43 1.43]                    ✅ FIB
T(8,5)    13    0.015103        0.27       [1.   1.59 1.59 2.08]                    ✅ FIB
T(5,3)    8     0.057954        0.50       [1.   1.   2.08 2.08]                    ✅ FIB

[2] WNIOSKI Z PORÓWNANIA
--------------------------------------------------------------------------------
Najbardziej harmoniczne węzły (TOP 3):
1. T(17,7) (Q=24) - MSE: 0.001914 
2. T(19,5) (Q=24) - MSE: 0.002056 
3. T(13,8) (Q=21) - MSE: 0.003680 ✅ FIB

[3] INTERPRETACJA FIZYCZNA
--------------------------------------------------------------------------------

Co to oznacza?
MSE mierzy, jak bardzo węzeł zachowuje się jak idealna, jednowymiarowa struna.
Węzły mocno splecione (jak T(13,11)) są "sztywniejsze" i mają zaburzone widmo
z powodu silnych oddziaływań geometrycznych między pętlami.
Węzły "luźniejsze" i asymetryczne (jak T(21,3)) mogą zachowywać się bardziej
jak swobodne pętle, co daje czystsze widmo rezonansowe.

Dla stabilnej cząstki potrzebujemy CZYSTYCH stanów kwantowych (długi czas życia).
Stany chaotyczne (wysokie MSE) szybko ulegają dekoherencji.

```
