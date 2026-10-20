# Nowy test teorii: niezależność statyczna nie zapewnia autonomii dynamiki

Data: 2026-10-20. Warunkowe badanie tej samej hipotezy
FIN-HYP-HIST-ME-001, bez dopisywania nowego postulatu fizycznego.

## Najważniejszy wynik

**Prawo wybrane dla większej całości na ogół nie jest prawem wybranym dla
samej obserwowanej części — nawet gdy miara stanów rozkłada się na iloczyn.**
To nowy test własnej reguły składania hipotezy, nie kolejna próba dopasowania
jej do HB. Odróżniono pominięcie współrzędnej istniejącego źródła od
zbudowania mniejszego modelu tą samą zasadą.

Po pominięciu L jednostajnych współrzędnych odniesienie dla pozostałych
zmienia się dokładnie z K na `K+L(q-1)I`. Niewidoczne zmiany nie znikają:
pozostają krokami własnymi w zredukowanym opisie. Nie dodano nowych pustych
znaków ani nie zmieniono mikroskopowego postulatu.

Dla q>=3 wyprowadziłem, że prawo następnego rzeczywiście innego stanu
pozostaje takie jak w izolowanym modelu **dokładnie przy jednostajnej
mierze obserwowanej części**. Dla niejednostajnej miary różni się również
wybór celu skoku, nie tylko zegar. W wymiernym przykładzie jest to
`1/3` zamiast `2/7`, z dokładną różnicą `1/21`.

## Drugi wynik: produkt równowagowy może dawać pamięć

Dla dwóch kopii o zadanej mierze `p_i nu_a` obserwowana pierwsza kopia
ma w tej hipotezie stacjonarną historię Markowa wtedy i tylko wtedy, gdy
ukryta miara nu jest jednostajna. Przy niejednostajnej nu powstaje pamięć,
choć w rozkładzie stanów nie ma korelacji między kopiami.

W dokładnym przykładzie dziewięciu stanów defekt dwukrokowego składania
wynosi `81/2560>0`. Jednocześnie kolejność rzeczywistych zmian pozostaje
jednostajna: w tym przykładzie pamięć dotyczy oczekiwania, a nie wyboru
następnej etykiety. Nie utożsamiam zatem każdej pamięci z niemożliwością
zmiany zegara.

## Zachowany przypadek dodatni

Przy jednostajnych miarach produktowych rodzina jest zgodna po jawnej
normalizacji pustych kroków. Wyprowadzono dokładne prawo pomijania kopii
i tożsamość próbkowania liczności bez zwracania. Podzbiór kopii jest w
dowodzie wybrany raz, a nie losowany od nowa przy każdym odczycie.

## Znaczenie i granice

Przechodzi **35 ograniczonych kontroli wymiernych**: źródła 9- i 27-stanowe,
zabite włókna, wariancja, wyjątek dwustanowy i przejście 10->6 typów
liczności. Kontrole te nie stanowią 35 dowodów ogólnych ani pomiarów.

Miary iloczynowe są kontrolami ogólnej hipotezy, **nie rzeczywistą
oddziałującą miarą A7 przy dodatnim gainie**. Nie wprowadzono niezależnych
części jako fundamentu nadsolitonu. Wynik mówi: jeśli od teorii wymagać
autonomii takich formalnie odłączonych sektorów, ta hipoteza na ogół jej
nie zapewnia. To wymaganie nie zostało przyjęte jako nowy aksjomat natury.

Pełny dowód: [PROOF_OR_OPEN_LEMMA.md](PROOF_OR_OPEN_LEMMA.md).
Nowy konkretny obiekt do dalszej analizy to rodzina `K+cI`, pochodząca z
już zadanego źródła. [Następny problem](NEXT_LIMIT_PROOF_ORDER.md) wymaga
osobnego dowodu granicy skalowanej; nie ogłaszam jej tutaj jako twierdzenia.
Strict F01 pozostaje otwarte, źródła i AGENTS są zachowane, cel badań aktywny.
