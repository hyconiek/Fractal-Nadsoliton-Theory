# ID-001: przyjęty budżet źródła i silniejsza granica błędu

Data: 2026-10-19. Wynik analityczny warunkowy, z ograniczonymi kontrolami;
nie formalizacja, niezależna recenzja naukowa ani walidacja fizyczna.

## Werdykt

**Przyjęto ograniczenie źródła o skończonej entropii dynamicznej dla
rzeczywistego, zadeklarowanego procesu N2. Nie wyprowadzono źródła dynamiki FIN.**

Źródłowy dowód poprawnie łączy Fano, stacjonarne sprzężenie całych historii,
regułę potęg entropii KS i dokładne prawo resetu N2. Zachowuje nierówne
prawdopodobieństwa p(g), pełną przekątną macierzy i zegar dwóch prób odświeżenia.
Współczynnik rbar=2(1-sum p²) jest częstością rzeczywistych zmian etykiety,
nie samą liczbą prób2. Entropie są w natach, nie w bitach starego ST383.

## Co nowego wyprowadziła recenzja

Nie poprzestano na sprawdzeniu liminf m e>=rbar. To samo sprzężenie można
obserwować co k kroków; jego błąd w chwili0 nie zmienia się. Otrzymujemy
nierówność F_12(e_(1/m))>=h_(k/m)-(k/m)H dla wszystkich k. Stąd przy jednym
skończonym H zachodzi **dodatnia granica samego błędu e**, nie tylko m e:

    liminf e_(1/m)>=F_12^-1(sup_(tau>0)[h_tau-tau H]_+)>0.

To mocniejszy wynik w tej samej klasie założeń. Nie ogranicza przybliżania
pojedynczych okien w słabszej metryce i nie mówi, że nadsoliton ma skończoną
entropię lub zewnętrzny czas. Dowód i kontrola g0 są w
[SUBSAMPLING_ERROR_FLOOR.md](SUBSAMPLING_ERROR_FLOOR.md).

Dodano również jawny błąd rozwinięcia entropii i doprecyzowano zastosowanie
do kompaktowego nośnika w R26. Sam skończony wymiar, histogram widma lub
nazwa chaosu nie zastępują miary, prawa ewolucji i właściwej entropii.

## Weryfikacja i pochodzenie

Źródło:14/14 manifestu. Z8 oryginalnych odwołań2 payloady zgadzają się
bajtowo,3 obecne odpowiedniki mają inną treść i3 oryginałów brakuje.
Nie zastąpiono ich fikcyjnymi zgodnościami. Bieżąca strategia47 plików
i87 wskazanych źródeł przechodzi kontrolę integralności w rzeczywistym repo.

Własny weryfikator ma23 ograniczone kontrole. Odtworzono źródłowy przykład
g3 na kopii; wyniki zgadzają się numerycznie, bez deklaracji identyczności
bajtowej. Maksymalny dryf zapisano w verification.json. Bez losowania
trajektorii, szerokiego skanu lub liczenia KS z histogramu. Wykorzystano
również aktualny generator144 stanów oraz odmienny przykład wymierny.

## Co pozostaje otwarte

F02 uzyskuje domknięcie **tylko jako przeszkoda klasowa**. F01 — źródło
stanu, konkretnego prawa, aktywności i kompozycji — nadal jest OPEN_SOURCE.
Nieskończony wymiar lub osobliwość nie obchodzą twierdzenia automatycznie;
trzeba wskazać rzeczywiście niespełnioną przesłankę. Dodatnia entropia
rekordów nie jest produkcją entropii: natywny model jest odwracalny
w równowadze i nie dostarcza w ten sposób strzałki termodynamicznego czasu.

## Kolejny łańcuch

Zgodnie ze strategią po ID-001 należy wrócić do decyzji źródłowej, nie
uruchamiać powtórnie ID ani kolejnego inwentarza. Jeśli pojawi się jeden
rzeczywisty obiekt historii/ewolucji, związać jego pochodzenie, zegar i odczyt,
zastosować ograniczenie entropijne i dopiero wyprowadzić jedną nową restrykcję
aktywności lub prawa. Bez obiektu gałąź pozostaje otwarta. Plan tej warunkowej
decyzji jest w FIN_POST_ID001_SOURCE_DECISION_20261019.md.

RP004-A/B/C pozostają gotową, osobno wybieraną gałęzią konsekwencji zadanego
działania. Ich ewentualny sukces nie zastąpi niezależnego źródła tego działania.
Nie dopisano nowego prawa, przestrzeni, zegara, kąpieli, jednostek lub
niższej warstwy informacji; brak promocji fizycznej, QM/GR/SM/ToE.
