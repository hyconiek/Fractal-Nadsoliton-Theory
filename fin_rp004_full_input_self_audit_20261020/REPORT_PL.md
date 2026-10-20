# RP004: audyt na pełnych lokalnych dowodach i jawne stałe

Data: 2026-10-20. Wykonano samodzielnie, zgodnie z odmową delegowania.
To audyt tego samego autora, nie niezależna recenzja ani formalna weryfikacja.

## Wynik matematyczny

Po przeczytaniu rzeczywistych lokalnych dowodów nie znalazłem błędu, który
unieważniałby warunkowy strumień RP004 w zamrożonym modelu. Podtrzymuję

    F_exact(R,T)>=(625/12386304) epsilon^4 T>0,
    T=12log(1/epsilon), R=2log(1/epsilon),

dla wszystkich dostatecznie małych epsilon, jednostajnie względem podanych
skończonych ograniczeń współczynników. R jest stałe w każdym całkowaniu
po czasie. To wynik dla zadanego klasycznego modelu 1+1, nie dowód emisji
fotonów, źródła działania, fizycznej czasoprzestrzeni ani ToE.

## Co sprawdziłem i doprecyzowałem

- Znak opóźnionej odpowiedzi, oba bieguny i normalizację Darboux; osobne
  sprawdzenie przez transformację Laplace'a i skok pochodnej Greena.
- Jednostajność stożka i pochodnych cząstkowych, razem z DC, modem kształtu,
  progiem i wysokimi częstościami. Jawne odcięcie 1/16 daje margines 1/12.
- H2 x H1 i ślady: wyprowadziłem bezpieczną stałą liniową 21 oraz konkretne
  stałe źródła i domknięcia bootstrapu. Nie opieram tego na H1 x L2.
- Wszystkie składniki oryginalnego prądu, oba końce przedziału i jeden próg
  epsilon_* działający dla całego zakresu mniejszych epsilon, nie tylko
  wybranego ciągu.

[Audyt przyczynowy](CAUSAL_AND_SPECTRAL_AUDIT.md),
[jawne stałe](EXPLICIT_STRONG_CONSTANTS.md),
[prąd i kwantyfikatory](SIGNED_CURRENT_AND_QUANTIFIERS.md).
Są to uzupełnienia sprawdzenia istniejącego twierdzenia, nie kolejny
niezależny wynik ani ponowne odkrycie dodatniego Gamma.

## Odtworzenie i pochodzenie

Nowy weryfikator przechodzi **47 ograniczonych kontroli**, wiąże 22 wejścia
i zachowuje 8 manifestów/115 wpisów poprzedników. Dwa stare skrypty naukowe
odtworzono tylko na kopii tymczasowej. Stary checker pakietu przechodzi
w zweryfikowanym kontekście historycznego prefiksu AGENTS, nie jako kontrola
dzisiejszego pełnego AGENTS. Osiem punktów kwadratury nadal nie jest dowodem
jednostajnego oszacowania stożka. [Wyniki](RESULTS.json), [replay](REPLAY.md).

Źródło zachowuje **10 zgodnych wpisów manifestu głównego i brak finalize.py**;
dwa manifesty zadań mają po 10 zgodnych wpisów. Wszystkie 33 obecne pliki
źródła pozostały bez zmian. W obecnym kontekście 7/10 oryginalnych odwołań
zgadza się w całości; trzy mają inne bieżące odpowiedniki. Jedno z nich,
AGENTS, ma dokładnie zgodny historyczny prefiks. To nie 10/10 oryginalnych
wejść ani dowód dawnego dostępu do późniejszych uzupełnień.

Stary audyt z brakującymi materiałami pozostaje niezmieniony. Dzisiaj pliki
były dostępne i przeczytane; nie naprawia to jego historycznej proweniencji.
AGENTS i dawne plany również zachowano.

## Co dalej

Ten etap kontroli RP004 jest zamknięty w zakresie samodzielnego audytu.
Niezależna recenzja nie została wykonana. Nie uruchamiam automatycznie
absorpcji, rozpraszania ani większych symulacji. **F01 pozostaje otwarte:**
poprawność konsekwencji zadanego działania nie wyprowadza jego źródła.
[Przekazanie i granice](HANDOFF.md).
