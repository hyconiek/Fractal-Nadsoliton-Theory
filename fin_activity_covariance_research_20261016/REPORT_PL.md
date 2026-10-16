# Wykonane badania: od swobody aktywności do odtwarzania szumu i nowego zadania

Data:2026-10-16. Status: badania analityczne z ograniczonymi kontrolami,
nie formalizacja ani niezależna recenzja naukowa.

Nie zatrzymałem pracy na ponowieniu OPEN_INPUT. W ramach obecnego zlecenia
wyprowadziłem trzy konkretne obiekty teoretyczne, a następnie zamroziłem
następne zadanie, które można od razu przekazać innemu agentowi.

## 1. Co wykonano

**Dokładne odtwarzanie aktywności z chwilowej kowariancji.** Dla znanej
miary Gibbsa liczności i jednokopiowych przejść pełny stanowo zależny tensor
Gamma wyznacza cały odwracalny generator. Odtwarzanie jest rekurencją na
włóknach liczności; ma jawne testy nieujemności i warunku końcowego. Nie
wymaga funkcji podziału. To silniejsze dane niż średnia kowariancja albo
jedna macierz szumu w punkcie krytycznym.

**Kontrprzykład silniejszy niż poprzednie przerzedzenie HB.** Zbudowałem
perturbację przewodności na sześciokątnym kole przejść, która zachowuje
Gibbsa oraz CAŁY punktowy dryf liczności, a zmienia aktywności, wariancje
chwilowe i czasy przebywania. Symetryzacja D12, dodatniość i jednostkowy
zegar prób mają jawne uzasadnienie. Zatem nawet dopisanie znajomości dryfu
do statyki nie wyprowadza jednoznacznej kinetyki.

**Tożsamość ustalająca średni szum i lemat krytyczny.** Jednakowa miara
stacjonarna i jednakowy dryf danego obserwabla wymuszają jednakową średnią
stacjonarną Gamma. Przy koncentracji szumu w lokalnym stanie FIN ustala to
współczynnik2 w krytycznym nawiasie martyngału. Nie dodano szumu ani prawa
ruchliwości; jest to skutek generatora, stacjonarności i jawnych warunków.

Dowody i ich dokładne granice: [THEOREMS.md](THEOREMS.md).

## 2. Weryfikacja i delegacja

Trzech agentów pracowało równolegle tylko odczytowo: nad perturbacją
zachowującą dryf, odwrotnością Gamma->aktywności i konsekwencją krytyczną.
Ich sprawdzenie pomogło doprecyzować warunki, lecz nie zastępuje dowodu.
Główny wykonawca napisał i uruchomił wspólny [weryfikator](verify_research.py).
Sprawdza rachunkiem wymiernym mały niejednorodny przykład q3/N3, q12/N3,
rekonstrukcję, odbite włókna i negatywne kontrole. Osobny przykład A7/g3
jest tylko diagnostyką binary64. Wyniki są w [RESULTS.json](RESULTS.json).
Małe N nie jest dowodem przebywania w asymptotycznym obszarze MT.

API rekonstrukcji odrzuca m.in. niecałkowite liczności, N0, puste dziedziny,
nieścisłe dane numeryczne i niedozwolone krawędzie. Nie jest to walidator
aparatury ani zamiennik wcześniejszych bramek przyjęcia dowodów.

## 3. Gdzie kończy się wykonanie

Brakujące niezależne źródło aktywności nie zostało wyprowadzone. Tensor
Gamma policzony z założonego HB nie jest niezależnym źródłem samego HB.
AS-002/003 i most fizyczny pozostają nieaktywne. Jest to wynik próby
wyprowadzenia, a nie zignorowanie zlecenia lub nieobecnego pliku.

Powstało natomiast gotowe nowe zadanie czysto teoretyczne: wykazać,
kiedy różne odwracalne kinetyki mają TEN SAM stacjonarny krytyczny proces
parzystości mimo swobody mikroskopowych aktywności. Pełnego twierdzenia
o tej klasie procesów nie ogłaszam jako już wykonanego.

Instrukcja wykonawcza:
[FIN_STATIONARY_KINETIC_RIGIDITY_RESEARCH_ORDER_20261016.md](../FIN_STATIONARY_KINETIC_RIGIDITY_RESEARCH_ORDER_20261016.md).
Agent ma udowodnić granicę procesu przy wystarczającym warunku koncentracji
szumu, a następnie sprawdzić, czy da się ten warunek osłabić. Ma oddać
dowód, poprawny kontrmodel albo precyzyjny brakujący lemat, nie kolejne wykresy.

Potencjał jest konkretny: mniejsza zależność teorii kolektywnej od arbitralnego
wyboru mikroskopowych aktywności. Nadal nie jest to źródło naturalnego A7,
przestrzeń, fizyczny zegar, energia zachowana, QM, grawitacja ani ToE.
Poprzednie pakiety i AGENTS pozostawiono bez zmian; nowy pakiet zawiera
tylko proponowany dopisek do przyjęcia po osobnym przeglądzie.
