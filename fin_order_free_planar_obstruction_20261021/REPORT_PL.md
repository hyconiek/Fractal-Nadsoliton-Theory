# Operator FIN: przeszkoda planarna niezależna od kolejności portów

Data: 2026-10-21. Nowy warunkowy dowód matematyczny i dokładny certyfikat;
nie niezależna recenzja, formalizacja ani wynik eksperymentalny.

## Wynik

Dotychczasowy certyfikat odrzucał bierną realizację kołowo-planarną tylko
w naturalnym porządku portów. Usunąłem tę zależność: **zadany pełny operator
L12 nie jest odpowiedzią dodatniej sieci rezystorowej w dysku dla żadnego
porządku jego dwunastu portów na pojedynczym brzegu**.

Wystarcza pięć portów: 0,1,2,3,6. Trzy dokładne nierówności wymuszają
konkretne pary przecinających się cięciw. Pierwsza ustala porządek 0,1,2,3;
dwie pozostałe umieszczają port 6 na rozłącznych łukach. To sprzeczność,
nie brak znalezionego porządku w numerycznej optymalizacji.

Dowód nie enumeruje 11! porządków. Weryfikator dodatkowo sprawdza wszystkie
12 porządków wybranych pięciu portów i dla każdego podaje ujemny minor.
Wagi cosinusowe i potęgi 9/5 mają nowe przedziały policzone wyłącznie
rachunkiem wymiernym, z jawnymi resztami Taylora i pierwiastkami całkowitymi.

## Odporność i znaczenie

Sprzeczność pozostaje prawdziwa, jeśli każda niediagonalna pozycja macierzy
odpowiedzi zmieni się najwyżej o 1/100000 w zadanej normalizacji. Zaburzenie
nie musi być cyrkulantą ani zachowywać D12. Nie jest to fizyczny budżet błędu;
pokazuje, że wyniku nie naprawia drobna korekta macierzy lub przestawienie
etykiet portów.

Nie istnieje też ciąg takich skończonych odpowiedzi zbieżny do badanego M,
nawet ze zmieniającą się kolejnością portów. Nie przenoszę tego automatycznie
na dowolne kontinuum: ewentualny most musi najpierw dowieść, że jego odpowiedź
rzeczywiście jest granicą sieci z tej klasy.

To ograniczenie konkretnej drogi wyłaniania lokalnej geometrii. **Nie jest
to zakaz wszelkiej lokalności FIN**: nie wyklucza sieci nieplanarnych,
wielu brzegów, mediów wektorowych/dynamicznych ani innych operatorów.
Badany L12 nie jest projekcją A7. Nie wyprowadzono fizycznej przestrzeni,
jej wymiaru, jednostek, kinetyki, QM, GR lub ToE; F01 pozostaje otwarte.

## Źródła i odtworzenie

Przez browser-harness przeczytałem w pierwotnym PDF definicje i warunek
minorów kołowych, z dokładnym sprawdzeniem konwencji znaków. Nagrywanie
pozostało wyłączone. Metadane i rzeczywisty hash PDF zapisano w
[LITERATURE_REGISTER.json](LITERATURE_REGISTER.json); pełny PDF nie jest
kopiowany do pakietu. Nie użyto Firecrawl ani delegacji.

[Pełny dowód](PROOF_OR_OPEN_LEMMA.md), [certyfikat](CERTIFICATE.json),
[wyniki kontroli](RESULTS.json), [replay](REPLAY.md).
Poprzednie źródła, manifesty, plany, AGENTS i istniejące zmiany zachowano.
Weryfikator ma 25 ograniczonych kontroli, wiąże 12 wejść i zachowuje 32
wpisy dwóch manifestów poprzedników. Te liczby nie oznaczają niezależnych
dowodów ani przeglądu wszystkich badań repozytorium.

Następna decyzja fundamentalna powinna uwzględnić zamknięcie całej tej
klasy pasywnego pojedynczego dysku, a nie szukać kolejnego porządku portów.
Inna klasa kompozycji lub lokalnego nośnika wymaga rzeczywistego, jawnego
obiektu i własnego testu; nie zostaje przyjęta automatycznie po tym no-go.
