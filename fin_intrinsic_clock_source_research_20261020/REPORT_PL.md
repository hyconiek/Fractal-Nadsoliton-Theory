# F01: wyprowadzono rozdział zegara i prawa następnego zdarzenia

Data: 2026-10-20. Nowe warunkowe badania analityczne, samoaudyt i małe
kontrole wymierne. Nie jest to niezależna recenzja ani domknięcie ToE.

## Wynik próby źródłowej

**Braku F01 nie da się w sprawdzonej klasie usunąć wyłącznie zmianą zegara,
nawet zależnego od całej dostępnej historii.** Pozostaje swoboda wyboru
tego, do którego stanu prowadzi następne zdarzenie. To precyzyjne nowe
ograniczenie próby wyprowadzenia, a nie kolejny audyt pustego wejścia.

Nie przeniesiono zadania do RP004: dodatni strumień w zadanej PDE nie
wyprowadza źródła jej prawa. Nie powtórzono filtra, G->Q ani kampanii
klasyfikatorów. Wykorzystano już istniejące ES/CF/RS.

## Nowe obiekty

1. **Dokładny podział kosztu całej historii.** Dla natężenia zmian a,
   warunkowych praw celu następnego zdarzenia r oraz docelowych s,

       koszt KL = koszt dopasowania zegara + integral a KL(r||s).

   Wyprowadzono także minimum po wszystkich dopuszczonych przyczynowych
   skalarnych zegarach. Pozostała część nie znika po zmianie jednostki
   czasu. Minimum jest narzędziem analizy, nie zasadą dynamiki FIN.
2. **Dokładne prawo pierwszego rzeczywistego wyjścia w ES.** Zsumowano
   wszystkie nieobserwowane odświeżenia i trzy fazy bloków. Otrzymane
   prawo nie jest automatycznie macierzą kolejnych skoków: start w chwili
   stacjonarnej różni się od startu na zdarzeniu. Tej różnicy nie pominięto.
3. **Jawna przeszkoda niezależna od zegara.** Pierwszy następny stan
   różni się od rodzimego resetu. Żadna zachowująca kolejność, pełna zmiana
   parametryzacji od tej samej chwili początkowej tego nie usunie.
   Osobno opłacono skończonoczasowe ograniczenie dla przyczynowej klasy
   zegarów; nie rozciągnięto go na zegary antycypujące przyszły cel skoku.

W istniejącym przykładzie wymiernym trzy etykiety dają dokładnie
`J_first(1,0)=1883/2502`, zamiast `3/4`. Dla modelowego T=1 wyprowadzono
`TV>1/3024` wobec **każdego** przyczynowego przeskalowania docelowej
kinetyki oraz minimalny koszt `KL>1/4572288`. To twierdzenia o modelu,
nie pomiary lub gwarancje praktycznych zasobów.

Ogólny dowód stosuje się również do rzeczywistego dodatniego, niejednorodnego
rozkładu N2 FIN i symetrycznie zdefiniowanego h z ES. Nie zbudowano jednak
pełnego odpowiednika 144-stanowego rodzica lub całej rodziny liczności.

Wynika stąd ścisła, ograniczona odpowiedź negatywna: w badanej klasie
rozkład stacjonarny, **wszystkie** prawa dwuczasowe i chwilowe kowariancje
nie wybierają nawet prawa historii z dokładnością do zmiany zegara.
Nie jest to zakaz przyszłego wyprowadzenia z rzeczywiście bogatszej,
niezależnie określonej organizacji procesów FIN.

## Co sprawdzono

Weryfikator ma **332 ograniczone kontrole**, w tym symboliczne logarytmy
liczb wymiernych, niezależną 81-stanową realizację ES, rozwiązywanie
27-stanowych zabitych włókien oraz dokładne marginesy 1728 słów q12.
Nie są to 332 dowody ani pełna enumeracja stanów FIN. Mały przykład nie
udaje rzeczywistego obszaru MT. Założenia nieskończonych szeregów i znak
wszystkich składników opłacono analitycznie, nie przez obcięcie numeryczne.

Zakres wyszukiwania wcześniejszych odpowiedników to 173 jawnie wybrane
pliki z istniejącego indeksu. DUPLICATE_SEARCH.json zachowuje zapytanie,
ścieżki i hashe. Brak trafienia w tym zapytaniu nie jest dowodem globalnej
nowości. Metody kompensatorów i entropii są odziedziczone i wskazane.

## Gdzie rzeczywiście zatrzymuje się wyprowadzenie F01

**Nie wyprowadzono fundamentalnego prawa ewolucji.** Wykazano natomiast,
że testowane dotychczas dane nie ustalają nawet jego części niezależnej
od zegara. Istnienie zmiany i korelacji nadsolitonu nie jest jeszcze
matematyczną regułą przypisującą warunkowe wagi przyszłym zdarzeniom.

Następny brak jest w [NEXT_OPEN_ATOM.md](NEXT_OPEN_ATOM.md). Jego treści
nie wolno przyjąć jako nowego aksjomatu tylko dlatego, że pozwoliłaby
otrzymać HB. CSP pozostaje nieaktywne; źródła, AGENTS i stare plany
pozostawiono bez zmian. Nie wykonano eksperymentu, delegacji ani commita.
