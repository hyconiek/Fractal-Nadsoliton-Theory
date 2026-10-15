# FIN SRN-001–SRN-005 — podsumowanie (tylko nowe badania)

Nowa obserwacja: przy każdej migawce z tej samej, z góry wybranej grupy liczymy dwie cechy (tryb 1 i 3) i zapisujemy wyłącznie, czy bezwzględna wartość trybu 3 jest większa niż trybu 1. Obie wartości można pomnożyć przez dowolną wspólną niezerową liczbę, nawet inną w każdym pomiarze: decyzja nie zmieni się. Nie działa to dla odrębnych wzmocnień, przesunięć ani dodatkowego szumu.

Nowy wynik analityczny: prawdopodobieństwa tego zdarzenia różnią się w dwóch modelach o d_ang*alpha + O(alpha²), gdzie d_ang≈0,02640349962 i alpha=k/N. Dowód używa rzeczywistych losowań bez zwracania elementów (hipergeometrycznych) i płaci prawdopodobieństwo leżenia blisko granicy |W3|=|W1|. W drugiej części użyto czterowymiarowego korektora Poissona i przeniesiono tylko końcowe ograniczone zdarzenie decyzji do niezmienionego FIN. Dla k=floor(N^(63/64)), M=floor(N^(1/4)) oba błędy klasyfikacji dążą do zera **warunkowo** względem istniejących przesłanek MT/KO/DKO/NOC.

Ważne ograniczenia: nie jest to prawo optymalnego kosztu, fizyczne urządzenie, eksperyment ani niezależny dowód doboru A7 przez przyrodę. Test nie jest odporny na osobne przeskalowanie kanałów. Wspólny mnożnik zera (zanik sygnału) również nie jest dozwolony.

**Uwaga praktyczna:** sekwencja jest skrajnie zachowawcza: warunek `k_N<=N/2` jest zapewniony dopiero dla astronomicznie dużych `N` (około `2^64`), a składnik `N^-1/64` maleje bardzo wolno. Nie należy z tych wykładników wyciągać wniosków o wykonalnym eksperymencie.
