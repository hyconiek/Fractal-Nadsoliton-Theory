# Ciągłe fazy nie wystarczają: ograniczenie liniowego źródła z mapą softmax

Data: 2026-10-21. Nowy warunkowy dowód dotyczący klasy źródeł. Bez nowego
prawa FIN, bez założenia skończonej dokładności i bez niezależnej recenzji.

## Co wyprowadziłem

Poprzedni wynik BI nie wykluczał bogatszej, ciągłej informacji początkowej.
Sprawdziłem więc konkretną klasę: ustalona liczba ciągłych współrzędnych,
liniowa/afiniczna propagacja i istniejąca w repo kanoniczna mapa
p(theta)=softmax(X7 theta), z odczytem znaku parzystości.

Nawet jeśli stan początkowy ma dowolną dokładność, dowolny rozkład i nośnik
fraktalny, liczba różnych binarnych historii dostępnych takim odczytem
rośnie tylko wielomianowo z długością zapisu. To twierdzenie o dostępie
do informacji, nie założenie, że ciągły stan zawiera jej mało.

Dokładniej, dla B odczytów i d ustalonych ciągłych współrzędnych jest
co najwyżej C(B+1)^d słów. Dopuszczenie dodatkowo12^N początkowych stanów
etykiet daje dostępną entropię najwyżej N log12+log C+d log(B+1).
Stała C nie zależy od wartości współczynników, z góry ustalonych chwil
próbkowania ani N;
nie obliczono jednak jej praktycznej wartości lub progu N0.

## Konsekwencja dla rzeczywistego modelu FIN

Rodzimy niepostselekcjonowany HB po opłaconym przygotowaniu ma na N^2
przedziałach odczytu (N^2+1 bitach) entropię rzędu N^2. Połączenie nowej
granicy pojemności z opłaconą nierównością BI/Fano daje dodatni próg
średniego błędu. **Samo ukrycie dowolnie dokładnej informacji w ciągłych
fazach nie odtwarza tych statystyk przez liniową propagację i kanoniczny
odczyt softmax.** Uwzględniono wspólne rozkłady etykiet i pola, bez ich
niezależności i bez stacjonarności konkurencyjnego źródła.

Wniosek obejmuje także wskazany analogiczny odczyt progowanego natężenia
klasycznych amplitud. Nie jest to zakaz liniowej mechaniki kwantowej z
losowym pomiarem Borna: taki pomiar nie jest deterministycznym progiem.

## Co sprawdziłem w literaturze

Przez browser-harness przeczytałem właściwe twierdzenia Wilkiego oraz
Aschenbrennera i współautorów. Kluczowe było zamienienie ról: współczynniki
odczytu są punktami zapytań, a stan początkowy jest parametrem rodziny.
Nowym wynikiem jest zastosowanie tych znanych narzędzi do rzeczywistej
mapy FIN i rodzimego budżetu zapisu, nie nowe twierdzenie o wymiarze VC.
PDF-y pozostają tylko w /tmp; metadane i rzeczywiste hashe są związane.

## Granice i dalsza luka

Nie wykluczam pełnej nieliniowej akcji, nieskończenie wielu losowych modów,
rosnącej złożoności odczytu ani rzeczywiście adaptacyjnego sprzężenia.
Kontrole przez mapę podwajającą i adaptacyjne progi pokazują, dlaczego nie
wolno rozszerzyć wyniku na wszystkie układy deterministyczne. Samo słowo
„nieliniowy” również nie dowodzi jeszcze właściwego prawa statystycznego.

Potrzebny obiekt źródłowy musi więc pokazać rzeczywistą złożoność dostępnych
historii i ich miarę, a nie jedynie powołać się na nieskończoną dokładność
stanu początkowego. Nie wyprowadzono F01, fizycznego zegara, detektora,
geometrii, QM/GR ani ToE. Nie wprowadzono niższej warstwy informacji pod
nadsolitonem. AGENTS i wcześniejsze źródła zachowano.

[Dowód](PROOF_OR_OPEN_LEMMA.md), [rejestr źródeł](LITERATURE_REGISTER.json),
[kontrole](RESULTS.json), [odtworzenie](REPLAY.md).
Przechodzą 23 nowe ograniczone kontrole LS oraz osobno odtworzenie 36
wcześniejszych kontroli BI. Związano 12 wejść i zachowano 40 wpisów
manifestów. Nie jest to komputerowy dowód twierdzeń Wilkiego/VC ani
niezależna recenzja analitycznej aplikacji.
