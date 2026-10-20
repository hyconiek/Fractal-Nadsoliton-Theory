# RP004: przegląd źródła, brakujący dowód i decyzja po strumieniu

Data: 2026-10-19. Nowe warunkowe dowody analityczne; nie formalizacja ani
niezależna recenzja naukowa. Bez zmiany źródeł, AGENTS lub praw FIN.

## Wynik

Uwzględniłem aktualny audyt strategii, `research_tasks`, rejestr blokad,
plan wykonania oraz późniejszy wynik ID-001. Źródło dynamiki pozostaje
najważniejszą otwartą zależnością. Wybrana kontynuacja RP004 rozwiązuje
inne pytanie: konsekwencję już zadanego klasycznego działania.

Nie skończyłem na przepisaniu OPEN. Wyprowadziłem brakujący obiekt:
**jednostajne oszacowanie odpowiedzi opóźnionej i pierwszych pochodnych
na rosnącym pasie T/6<|x|<T/3, z błędem O(T^-1/2)**. Dowód oddziela
biegun harmoniczny od przejściowych punktów stacjonarnych, zachowuje DC,
włączenie wymuszenia referencyjnego, mod kształtu i oba kierunki widmowe.

Naprawiłem też błąd domknięcia stałych w dostarczonym bootstrapie B.
Najpierw poprawia się oszacowanie r, dopiero potem w; nie można ukryć
początkowego współczynnika przy epsilon X w dodatkowej małości.

Łącznie daje to dla oryginalnego nieliniowego rozwiązania

    F_exact(R,T) >= (625/12386304) epsilon^4 T >0,
    T=12log(1/epsilon), R=2log(1/epsilon),

dla dostatecznie małego epsilon w zadanym zakresie a i skończonych granic
współczynników. Promień pozostaje stały podczas całkowania po czasie.
Nie dodano początkowego pakietu falowego, Phi6 ani zewnętrznego sterownika.
Nie jest to dowód tłumienia przez cały czas, fotonu lub fizycznego źródła FIN.

## Uczciwy stan materiałów

Źródło samo pozostawiło A-GATE otwarte i nie wykonało C. Nowy dowód jest
naszą pracą tutaj, nie przypisanym mu historycznym sukcesem. Z 31 pozycji
trzech manifestów 30 odpowiada obecnym plikom, a `finalize.py` brakuje.
Nie odbudowuję nieobecnego skryptu ani nie twierdzę, że go sprawdzono.
Osiem z dziesięciu unikatowych
referencji wejściowych pasuje; dwie mają inne lokalne odpowiedniki.

Oba skrypty źródła odtworzono tylko na kopii. Nowy weryfikator wykonuje
26 ograniczonych kontroli i osiem diagnostyk całki przyczynowej. To kontrola
algebry, znaków i możliwych błędów wzorów, **nie 26 dowodów PDE**. Wszystkie
33 obecne pliki źródła i historyczny snapshot strategii są zachowane.

## Następny krok

Gotowa [decyzja i instrukcja](../FIN_POST_RP004_RESEARCH_DECISION_20261019.md)
zleca najpierw jeden niezależny audyt nowego dowodu i budżetu reszt.
Nie zleca ponownego liczenia Gamma, symulacji promieniowania ani absorpcji
na drugiej ścianie. Niezależnego audytu nie wykonano w tej sesji.

Dla kierunku fundamentalnego nadal potrzebny jest rzeczywisty, niekołowy
obiekt ustalający prawo historii, aktywność, nieliniowość lub kompozycję.
Bez niego zachowujemy OPEN_SOURCE, a nie nowy pusty inwentarz. Modelowy
transport energii nie ustala źródła A7, zegara, geometrii, jednostek, QM/GR
ani ToE. Aktywny przykład nadal nie bada bloku spektatorów K23.
