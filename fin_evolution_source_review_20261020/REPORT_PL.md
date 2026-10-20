# ES wcielone; dalsze wyprowadzenia: filtr, całe trajektorie i defekt pamięci

Data: 2026-10-20. Badania warunkowe, bez zmiany źródeł i bez nowego prawa FIN.

## Co wcieliłem

Przyjmuję poprawne wnioski ES: zgodność skończonej hierarchii korelacji oraz
Gamma nie wymusza zgodności całej historii bez uzasadnionej Markowowskości.
Kryterium martyngałowe jest poprawne, lecz samo nie stanowi niezależnego
źródła ewolucji. **F01 nie zostało zamknięte.**

Źródło zachowano: 15/15 pozycji manifestu, wszystkie 16 plików i 14 wejść.
122 kontroli producenta odtworzono tylko na kopii; wynik jest tu bajtowo
zgodny. To odtwarzalność, nie zewnętrzna recenzja lub 122 dowody.

## Co wyprowadziłem dalej

1. **Dokładny stan predykcyjny.** Dla zadanego generatora źródła, odczytu
   i przygotowania wyprowadziłem filtr warunkowy. Uwzględnia on ukryte
   przejścia, które nie zmieniają obserwowanej etykiety. Nie zakładamy już,
   że sama bieżąca etykieta jest pełnym stanem.
2. **Błąd całego prawa trajektorii.** Całka z różnicy rzeczywistych szybkości
   przewidywanych z historii ogranicza TV całej ścieżki. Wyprowadziłem też
   dokładną tożsamość kosztu entropijnego, z zachowaniem błędu początkowego,
   i najlepsze stałe streszczenie Markowowskie. W ES jest nim właśnie B,
   lecz dodatni koszt pominiętej pamięci pozostaje.
3. **Jawny defekt pamięci.** Dla rodziny ES obliczyłem dodatnią dolną granicę
   tego defektu mimo identycznych korelacji par. Jednocześnie pokazałem,
   dlaczego nadmiarowa wiedza pełnego źródła może dawać niezerowy surowy
   defekt nawet przy dokładnie Markowowskim odczycie. To ważne rozróżnienie
   między warunkiem wystarczającym a koniecznym.
4. **Prawo oczekiwania.** W wymiernym przykładzie różnica praw braku zmiany
   zaczyna się od `t^2/96`, choć wszystkie dwuczasowe końce przedziału się
   zgadzają. Jest to wynik o kontrmodelu ES, nie zmierzona aktywność natury.

Dowody: [filtr](EXACT_CAUSAL_FILTER.md),
[transfer pełnej ścieżki](FULL_PATH_TRANSFER.md),
[kontrprzykłady i dolne granice](FILTERED_BLOCK_DEFECT.md).
Nowy weryfikator ma 83 ograniczone kontrole, głównie wymierne; nominalna
kontrola rodzimego N2 jest oddzielona od dowodu. Nie użyto nowego szumu,
kąpieli, sterownika ani fizycznego zegara do naprawienia modelu.

## Co pozostało naprawdę otwarte

Potrzebny jest niezależnie uzasadniony model źródła wraz z mapą odczytu
i przygotowaniem. Wtedy obecne twierdzenia dają już konkretny sposób
wyprowadzenia prawa odczytu i oszacowania błędu, zamiast dopisywania jego
Markowowskości. Takiego fundamentalnego źródła nadal nie dostarczono.

[Gotowy następny punkt badawczy](../FIN_CAUSAL_SOURCE_TRANSFER_NEXT_CHAIN_20261020.md)
wiąże dokładnie ten brak. Nie zleca kolejnego katalogu korelacji ani pustego
inwentarza. Źródła, AGENTS i poprzednie pakiety są zachowane; bez commita.
