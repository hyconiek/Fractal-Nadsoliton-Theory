# AH-001--AH-003: przyjęcie pełnej historii anonimowej i skończonego budżetu

Data: 2026-10-20. Integracja w drzewie roboczym, bez commita i bez zmiany
źródła, AGENTS lub wcześniejszych pakietów.

## Wynik

Przyjmuję dwa główne wyniki warunkowe:

1. **Pełne anonimowe historie identyfikują łańcuch do zmiany nazw stanów.**
   Dotyczy to zadanego dodatniego 12-stanowego modelu skokowego. Dowód
   zachowuje losowe nazewnictwo pierwszych odwiedzin i zależność przejść.
2. **Istnieje jawny wystarczający skończony prefiks.** Przy znanej wspólnej
   dolnej granicy a dla obu modeli, gwarantowanej separacji delta i zadanym
   błędzie eta, podany L zapewnia błąd decyzji najwyżej eta pod każdą hipotezą.

Doprecyzowałem kwantyfikatory i filtrację, samodzielnie wyprowadziłem
nierówność martyngałową używaną w oszacowaniu i sprawdziłem stałe. Dodatkowy
argument z istnienia skończonego testu niezależnie potwierdza implikację
identyfikowalności, zamiast polegać wyłącznie na granicy nieskończonego rekordu.

## Co rzeczywiście sprawdzono

Manifest **15/15**, zachowane **16 plików**. Sześć z siedmiu referencji
wejściowych zgadza się; oryginalnego wklejonego polecenia brakuje. Oba
skrypty źródła odtworzono wyłącznie na kopii; wynik JSON jest tu bajtowo
zgodny. Własny weryfikator ma **49 ograniczonych kontroli**, nie 49 dowodów.

Rachunek wymierny z kontrolą reszty szeregu logarytmu potwierdza przykładowe
`m=569231067`, `L=1138462136`. To ponad miliard **obserwowanych stanów**,
nie praktyczne zalecenie pomiarowe. Nie symulowano takiej historii. Parametry
tego wzoru są oddzielne od przykładu 3200 przejść z wagami 1,...,6.

## Następny krok

Anonimową rekonstrukcję zamykam w jej obecnym zakresie. Nie zalecam kolejnej
kampanii motywów lub poprawiania wykładnika. [Gotowy dalszy łańcuch](../FIN_POST_AH_SOURCE_AND_VALIDATION_NEXT_CHAIN_20261020.md)
rozdziela istniejący audyt RP004 od fundamentalnego przyjęcia rzeczywistego
nowego obiektu źródłowego. Bez takiego obiektu nie uruchamiać pustego inwentarza.

Wynik mówi, jak odtworzyć **założone** prawo z bogatej historii, a nie dlaczego
nadsoliton miałby takie prawo wybrać. F01/F03/F06, fizyczny zegar, geometria,
aparatura/EV oraz QM/GR/ToE pozostają otwarte. Model jednolabelowy nadal nie
jest rodzimym procesem A7/N2 ani teorią polową.
