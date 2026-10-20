# Oryginalne źródło hipotezy ma już opłacone okno lokalnej równowagi

Data: 2026-10-20. Dowód warunkowy i samodzielna kontrola, nie niezależna
recenzja ani wyprowadzenie prawa natury.

## Co zmieniło się naukowo

Bez dopisywania jednostajnych współrzędnych ukrytych udowodniłem
wielomianowe mieszanie lokalne oryginalnego źródła entropijnego na N
oddziałujących współrzędnych oraz jego okno lokalnej równowagi po starcie
iid. Zastosowano rzeczywistą globalną miarę Gibbsa A7, pełny warunek
skalowania i właściwe czynniki liczności. Nie wybrano nowego lokalnego
optymalizatora zamiast zadanego globalnego źródła.

Kluczowy nowy lemat wykorzystuje pełne kliki zmian jednej współrzędnej.
Daje globalny iloraz sąsiednich skalowań O(N), prawdopodobieństwo pozostania
co najmniej rzędu N^-2 i aktywność co najmniej rzędu N^-1. To wystarcza
do porównania form, także z opłaceniem ujemnego widma w czasie dyskretnym.

**Nie udowodniłem mocniejszej stałej dolnej granicy wymaganej w dawnym
zleceniu.** Zachowałem ten brak. Jawny model kontrolny pokazuje ponadto,
że same ograniczone ilorazy wag i symetrie nie wystarczają do takiej stałej.
Jego miara nie jest jednak Gibbsem A7, więc nie jest to kontrprzykład dla
właściwego FIN. Słabszego wyniku nie przemianowano na mocniejszy.

## Co nadal trzeba wyprowadzić

Granica krytyczna oryginalnego źródła pozostaje otwarta. Potrzebne są dwa
konkretne oszacowania jego rzeczywistego dryfu i surowej wariancji. W czasie
dyskretnym nawias martyngału wymaga jeszcze odjęcia kwadratu średniego
przyrostu; wyprowadziłem tę poprawkę i warunek, który ją kontroluje.
[Następne zlecenie](NEXT_GLOBAL_ACTIVITY_DEFECT_ORDER.md) wskazuje dokładne
wielkości, zamiast ponawiać pustą bramkę źródłową lub serię symulacji.

To wynik w upoważnionej hipotezie FIN-HYP-HIST-ME-001. Nie uzasadnia samej
zasady maksymalizacji entropii jako prawa natury. Naturalne F01/F03/F06,
źródło A7, fizyczny czas, jednostki, geometria, kalibracja, QM/GR i ToE
pozostają otwarte. Źródła, AGENTS i historyczne plany są zachowane.

Szczegóły: [dowody](PROOF_OR_OPEN_LEMMA.md), [kontrole](RESULTS.json).
Przechodzi 67 ograniczonych kontroli dokładnych, nie 67 dowodów ogólnych.
Związano 15 plików wejściowych i zachowano 292 wpisy trzynastu manifestów.
