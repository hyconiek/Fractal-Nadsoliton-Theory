# Przyjęcie pięciu badań po RP004 i korekta interfejsu kowariancji

Data: 2026-10-20. Integracja w drzewie roboczym, bez commita i bez zmiany
oryginałów, AGENTS lub historycznych planów.

## Co przyjmuję

Pięć pakietów daje poprawne twierdzenia **warunkowe**: dokładną pamięć po
eliminacji połowy zmiennych, ograniczenia działania KL/Fishera, rekonstrukcję
sześciu szybkości z oznaczonych korelacji oraz rozróżnianie praw przez historie
zdarzeń. To identyfikacja i konsekwencje założonego prawa, nie jego źródło.

Szczególnie konkretny jest świadek dwóch nieizomorficznych generatorów,
które mają takie same powroty dla wszystkich opóźnień. Cztery kolejne stany
rozróżniają je przez `ABAB`, z dokładną różnicą `11/34629528` i TV pełnych
partycji `11/17314764`. Cztery stany oznaczają tu **trzy skoki**. Nie wynika
stąd identyfikowalność dowolnego prawa z czterech stanów.

## Co poprawiam

Źródło zbyt ogólnie odrzuca uśrednioną kowariancję. Wyprowadziłem i sprawdziłem
rachunkiem wymiernym korektę: dla odwracalnego jednolabelowego procesu
one-hot i znanej dodatniej pi, pełna średnia macierz `G` już wyznacza

    Q=-(1/2)diag(pi)^(-1)G.

Nie trzeba wtedy znać wszystkich warunkowych `Gamma_i`. Nie przenosi to
wyniku na proces liczności FIN, nie odzyskuje prądów nieodwracalnych i nie
ustanawia niezależnego źródła G. [Dowód i granice](MEAN_BRACKET_CORRECTION.md).

## Co sprawdziłem

Manifest zbiorczy: **75/75**; zachowano wszystkie **77 obecnych plików**.
Odtworzono pięć skryptów naukowych tylko na kopii; własny weryfikator ma
**38 ograniczonych kontroli**, w tym dokładną enumerację partycji i rachunek
w Q(sqrt(3)). Różnice numeryczne sięgają około `3.82e-14`, nie są dowodem
błędu twierdzeń ani zgodnością bajtową. Brakujące oryginalne ZIP-y i dwa
wklejone dokumenty pozostają jawnie brakujące. Nie przypisuję checkerowi
nieobecnej w kodzie weryfikacji ZIP CRC.

## Następny łańcuch

[Gotowa instrukcja](../FIN_RELATIONAL_HISTORY_NEXT_CHAIN_20261020.md) rozdziela:

1. istniejący, wykonalny niezależny audyt RP004 — pliki dowodu są już lokalnie;
2. opcjonalne nowe badanie: czy **pełna** anonimowa historia identyfikuje
   łańcuch do przemianowania stanów i jaki skończony fragment wystarcza;
3. fundamentalny most: dopiero rzeczywisty niezależny obiekt źródłowy, nie
   kolejne liczenie historii wygenerowanych z założonego Q.

Audytu RP004 ani nowej kampanii anonimowych historii tutaj nie wykonano.
F01/F03/F06 oraz jednostki, geometria, walidacja fizyczna, QM/GR i ToE
pozostają otwarte. Nie utożsamiam modelu Z12 z oddziałującym A7/HB lub PDE.
