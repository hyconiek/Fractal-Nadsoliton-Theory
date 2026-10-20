# MH-001: opłacona jedna rosnąca rodzina przekrojów Halla

Wykonano: 2026-10-20. Nazwa katalogu zachowuje datę pierwotnego zlecenia.
Dowód analityczny warunkowy względem przyjętych A7/MT/UF/PC; ten sam autor,
bez formalizacji, subagentów i niezależnej recenzji naukowej.

## Nowy wynik

**ONE_MESOSCOPIC_FAMILY_PAID.** Dla zamrożonej kuli o promieniu 1/10 we
współrzędnych szybkich, wszystkich dostatecznie dużych N podzielnych przez
12 i sąsiednich poziomów w [-1,1] opłaciłem wymagane SHELL i ADJACENT.
Z dokładnej tożsamości dopełnień wynika warunek Halla dla usunięcia tej kuli
z całego lewego włókna. Jej liczba stanów rośnie jak N^5, a masa warunkowa
nie zanika. Nie jest to kolejna rodzina skończenie wielu dziur.

## Jak usunięto rzeczywistą trudność

- Wyprowadziłem dokładnie sześć przesunięć siatki A5+A5, jej kowolumen 6
  oraz rozkład na dwuwymiarowy rzut fazowy i ośmiowymiarowe jądro. Zostają
  skończone klasy reszt; nie zakładam ortogonalnego iloczynu całej siatki.
- Dla ostrej kuli ważone sumowanie Poissona daje błąd O(N^-10/11).
  Jest mniejszy niż szerokość potrzebnej powłoki N^-1/2. Jej dodatnia masa
  nie została wywnioskowana wyłącznie ze słabej granicy Gaussa.
- Normalizacja fazy wymaga osobnego dwuwymiarowego oszacowania dysku.
  Po dokładnym rozkładzie 8+2 jego błąd wynosi O(N^-2/3), też o(N^-1/2).
- Użyłem gładkiego przedłużenia rzeczywistych wag liczności, wraz z
  silnią/multiplikacją, oraz dokładnego znoszenia prądów i źródła między
  sąsiednimi poziomami. Prawdziwa warstwa wyjścia E_N i jej poprawki są
  jawnie opłacone, nie usunięte z modelu.

[Pełny dowód](PROOF_OR_OPEN_LEMMA.md), [SHELL](SHELL_PROOF_OR_OPEN_LEMMA.md),
[ADJACENT](ADJACENT_PROOF_OR_OPEN_LEMMA.md),
[dokładny zakres Halla](MESOSCOPIC_HALL_DISPOSITION.md).

## Weryfikacja i literatura

Przechodzą **33 ograniczone kontrole** dokładnego rachunku wymiernego:
siatki, 54 dozwolonych klas przesunięć jądra, sześciu fixture'ów liczności
i wszystkich 36 ruchów na fixture, wykładników błędów, znoszenia źródła
i poprawki E. To nie 33 dowody asymptotyczne ani enumeracja całego grafu.
Małe przykłady N120 badają algebrę, nie rzeczywistą dziedzinę MT.
Związano 17 wejść i zachowano 202 wpisy dziesięciu manifestów poprzedników.
[Wyniki](RESULTS.json), [odtworzenie](REPLAY.md).

Zgodnie z Twoją prośbą użyłem **browser-harness**, nie Firecrawl. Połączenie
z rzeczywistą stroną arXiv i pobranie PDF przez CDP są potwierdzone; nagrywanie
pozostaje wyłączone. Przeczytałem odpowiednie fragmenty pracy Guo o sumowaniu
Poissona i wygładzaniu. Ogólna metoda jest znana; nowe jest zastosowanie
do tej konkretnej siatki, wag i prądów FIN. [Zakres lektury](LITERATURE_NOTE.md).

## Granice wyniku i powrót do architektury

To nadal NIE jest pełny Hall dla dowolnych zbiorów, sprzężenie wszystkich
prądów, generator z cap5N, sklejenie całego D_N ani kontrprzykład do IBA.
Zachowuję twierdzenie KR pod NC/IBA i otwarty status bezwarunkowego IBA.
Nie uruchamiam kolejnych promieni lub wariantów kuli.

Wynik dotyczy odporności matematycznej zadanego rodzimego modelu HB,
nie niezależnego źródła kinetyki. F01, naturalne A7, jednostki, geometria,
walidacja fizyczna i QM/GR/ToE pozostają otwarte. Oryginały, AGENTS i
historyczne plany zachowano. [Handoff](HANDOFF.md).
