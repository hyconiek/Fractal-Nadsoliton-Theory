# Poszukiwanie źródła ewolucji: wynik i dokładnie odsłonięta przesłanka

Data: 2026-10-20. Nowe ograniczone badanie teoretyczne, nie recenzja zewnętrzna,
formalizacja lub przyjęcie nowego prawa FIN.

## Odpowiedź na zlecenie

**Nie odnalazłem niezależnie wyprowadzonego źródła ewolucji z obecnych
przesłanek FIN. Nie ogłaszam zamknięcia F01.** Nie oznacza to braku dynamiki
w modelach: generator HB i działanie polowe definiują rzeczywiste prawa
warunkowe. Nadal nie mamy ich pierwotnego wyboru z organizacji nadsolitonu.

Nie poprzestałem na ponownym spisie braków. Sprawdziłem mocniejszy kandydat
na obiekt źródłowy: pełne korelacje dwu-czasowe i chwilowe kowariancje
rodzimego procesu różnicy N2, bez dopisanej Markowowskości. Powstała jawna
rodzina diagnostycznych historii, które zachowują te dane **dokładnie**,
a mimo to mają inne prawo całej historii. Zjawisko ogólne jest znane w
literaturze; nie przypisuję FIN jego odkrycia.

## Co wyprowadzono

Z dodatniego rozkładu p konstruuję dodatnie bloki

    T_r=prod p_i [1+epsilon prod h_i], E_p h=0,

a następnie stacjonarną historię z losową fazą bloków i porównawczym
zegarem Poissona. To **kontrmodel analityczny**, nie nowa warstwa informacji,
kąpiel lub przyjęte prawo nadsolitonu.

Dla r=3 konstrukcja zachowuje wszystkie macierze

    P_t=exp(-2t)I+(1-exp(-2t))1p^T,

oraz pełne przyrostowe Gamma_i warunkowane tylko bieżącą etykietą i średnią G.
Ale jej trójczasowy moment różni się od resetu Markowa o

    (epsilon/3)(E_p h^2)^3 (2t)^2 exp(-4t).

Nie jest to tylko równość widma lub jednej funkcji powrotu. Dla każdego
ustalonego skończonego rzędu k można zachować całe k-czasowe prawa i zmienić
wyższe historie. Nie należy więc ratować bramki kolejnym skończonym katalogiem
korelacji ani uznawać zgodności macierzy P_t z semigrupą za dowód Markowowskości.

Wyprowadziłem też pozytywne kryterium: **pełnohistoryczne martyngały
współrzędnych** z autonomicznym B wymuszają rzeczywiste prawo Markowa
`exp(tB)`. Dowód przez równanie wsteczne jest w [PROOF.md](PROOF.md).
Nie wyprowadziłem jednak tych martyngałów z ontologii FIN. Wpisanie do nich
pożądanego B byłoby kołem, nie znalezieniem źródła.

## Co to zmienia w F01

Najbliższą potrzebną przesłanką nie jest następna liczba Gamma. Trzeba
wyprowadzić **jaki stan nadsolitonu jest kompletny dla przewidywania** oraz
czy jego rzeczywisty odczyt ma autonomiczne prawo warunkowane całą przeszłością.
Dopiero w opłaconej klasie można użyć istniejących odwrotności G->Q lub
historie->J do wykluczania alternatywnych aktywności. Skończona projekcja
pełnego stanu może mieć pamięć; nie wolno narzucić jej Markowowskości
dlatego, że chcemy odzyskać HB.

[Konkretny kontrakt brakującego atomu](SOURCE_ATOM.md) podaje wymagane mapy
i warunek do wyprowadzenia. Nie jest nowym aksjomatem ani uruchomieniem
kolejnego pustego inwentarza. Nowy postulat wymaga osobnego oznaczenia i
zgody autora; obecna praca go nie przyjmuje.

## Weryfikacja i granice

Weryfikator sprawdza dokładne dodatnie bloki, właściwe marginesy, przesunięcie
i odwrócenie skończonych okien, kontrast wyższych historii i kowariancje przyrostów
rachunkiem wymiernym. Osobno wiąże nominalny p z wersjonowanym safe_seed.
To ograniczone kontrole, nie mechaniczny dowód wszystkich twierdzeń ani
zewnętrzna walidacja. Pełne szczegóły i liczba kontroli są w RESULTS.json.

Źródła, AGENTS i dotychczasowe pakiety zachowano. F01/F03/F06, fizyczny
zegar i geometria, strict/legacy/selektor/role, B4/FIN-M0 oraz QM/GR/SM/ToE
pozostają otwarte. Nie wykonano audytu RP004 ani nowego eksperymentu.
