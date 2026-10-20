# Aktualizacja łańcucha ES/CF po przeszukaniu istniejących badań

Data: 2026-10-20. Nowa praca analityczna i małe kontrole, bez nowego
postulatu FIN i bez niezależnej recenzji naukowej.

Zamiast powtarzać wyprowadzenie filtra lub inwersję G->Q wykorzystano:
HM: dodatni defekt kompresji; CF: pełnohistoryczne intensywności i błąd
ścieżki; ES: pamięć niewidoczną w całym prawie dwuczasowym; P276: test
zgodności dwóch zadanych macierzy z jednym generatorem. To różne zakresy.
Zachowano też wcześniejszą bramkę bariery aktywacji: odpowiedni pośrednik
śledzony w zegarze śladowym jest konstrukcją zgodności, nie źródłem HB.

Nowy wynik: dla **zadanego odwracalnego źródła Markowa**
`delta=tr(C(2t)-C(t)^2)` jest dokładną nieujemną normą ukrytego defektu.
Daje `TV<=M sqrt(delta)/2` dla całego rekordu M przejść na tej siatce.
Przy dodatnich przejściach dostępne jest również oszacowanie KL.
Dokładne zero przy jednym dodatnim opóźnieniu wymusza dla skończonego
odwracalnego CTMC pełną zgodność Markowa, także poza siatką.

To istotnie zawęża interpretację ES: odwracalność obserwowanej historii
nie oznacza odwracalności jej ukrytego źródła. Tego założenia nie wolno
dopisywać do nadsolitonu. Źródłowa reprezentacja ES jest skierowana.

Zachowano negatywny wynik: mały defekt na jednej rzadkiej siatce nie
kontroluje całej ścieżki ciągłej. Jawny odwracalny przykład trójstanowy
ma taki defekt dążący do zera, a błąd ścieżki pozostaje większy niż 7/300.
To kontrola zakresu twierdzenia, nie nowy kontrmodel q12/IBA.

Pełne dowody: [PROOF_OR_OPEN_LEMMA.md](PROOF_OR_OPEN_LEMMA.md).
Zakres grepów i trafienia: SEARCH_RESULTS.json po wykonaniu poniższego
replayu. Nie twierdzimy, że przeczytano każde przeszukane źródło ani całe
repo. Nowość dotyczy tej lokalnej integracji narzędzi, nie pierwszeństwa
matematycznego w literaturze.

F01 nie została zamknięta. Następny most źródłowy wciąż wymaga rzeczywistego,
niekołowego dostawcy prawa/historii i mapy odczytu. Nie należy uruchamiać
CSP wyłącznie na niniejszej nierówności ani ponawiać pustego inwentarza.
Oddzielnie pozostaje audyt właściwych lokalnych dowodów RP004; ta praca
go nie wykonuje i nie promuje dodatniego strumienia do wyniku niezależnie
zweryfikowanego.
