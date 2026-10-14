# Podsumowanie badań TH-001–TH-005

Najmocniejszym nowym wynikiem jest pełne zastosowanie twierdzenia Kraaija o dużych odchyleniach do zadanej dynamiki FIN. Otrzymany koszt trajektorii nie jest przyjętą zasadą działania, lecz wynika z generatora skokowego i ma postać kosztu przepływów Poissona po przejściach etykiet.

Drugi ścisły wynik pokazuje, że siedem średnich modów nie wystarcza do przewidywania fluktuacji trajektorii: dwa prawdziwe stany licznikowe N=120 mają tę samą średnią, lecz różny generator kwadratu obserwabli. Cztery dodatkowe momenty drugiego rzędu wystarczają natomiast do dokładnego odtworzenia całego wektora zajętości p.

W obszarze folda wyprowadzono z rzeczywistego generatora współczynnik szumu miękkiego modu: jego wariancja na jednostkę modelowego czasu wynosi dokładnie `2 g_f/N`, więc po krytycznym skalowaniu amplituda graniczna to `sqrt(2 g_f)`. Nie jest to szum dopisany na poziomie makro.

Skalowanie `g-g_f~N^-2/3`, amplituda `N^-1/3` i czas `N^1/3` prowadzą formalnie do zatrzymanej dyfuzji

`dZ=-g_f(a lambda+b Z^2/2)ds+sqrt(2g_f)dB`.

Nie ogłaszam jednak pełnego twierdzenia zbieżności. Pozostał jeden konkretny lemat: trzeba rygorystycznie wyeliminować szybkie mody stabilne i ukryte składowe zajętości w problemie martyngałowym na przyspieszonym czasie. Diagnostyka skończonego N wspiera dokładnie przewidywane poprawki, ale nie zastępuje tego dowodu.

## Weryfikacja pakietu

Niezależny skrypt odbudował model bez importowania kodu producenta i przeszedł 41/41 kontroli. Czysty replay odtworzył kluczowe pliki wynikowe bajt w bajt. To kontroluje rachunki i implementację, ale nie zastępuje opublikowanego dowodu twierdzenia LDP ani brakującego lematu TH-005.
