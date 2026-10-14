# Wyniki badania pamięci FIN – 14 października 2026

Badanie dotyczy **modelu skończonego** (3 kopie, 12 etykiet, 1728 stanów), a nie empirycznego potwierdzenia fizyki FIN.

- D01: zgodne prawo stacjonarne i generator tego samego rodzica.
- D02: referencyjna propagacja Poissona, dane wejściowe w przedziałach `mpmath.iv`, modelowany budżet referencji 1e−8 TV (warunkowy numerycznie).
- D03: jawna pamięć `K(t)` i wymuszenie od nieekwilibrialnego przygotowania `b(0)`. Obliczono oddzielnie K(0), K'(0) i numeryczne momenty całkowe M0/M1. Uniwersalne zamknięcie ze stałą bazą NIE zostało udowodnione.
- D04: przed otwarciem odłożonych wyników zapisano algorytm, predykcje 144/288 elementów, oszacowania i ich SHA-256.
- D05: **36/36** testów dla FIN oraz dwóch kontrmodeli spełniło TV ≤0.01 i zadeklarowane oszacowania. FIN max TV modelu z pamięcią **2.733826296e-08**, a analogicznego modelu Markowa **0.3548376953**.

Ważne: dobra redukcja jest **warunkowana przygotowaniem i kanałem odczytu**; dla nowej preparacji baza Kryłowa jest generowana ponownie według niezmiennej reguły bez dopasowywania parametrów. Nie jest to odkrycie jednego niezależnego od stanu procesu Markowa ani dowód unikalności A7. Zgodnie z instrukcją RG-D06 zostaje zarezerwowany do decyzji architektonicznej.
