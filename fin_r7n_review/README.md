# Audyt i integracja pakietu R7N — 2026-09-20

Źródło: [FIN_R7N_HANDOFF_20260920](../FIN_R7N_HANDOFF_20260920/HANDOFF.md).
Oryginalny pakiet pozostaje niezmieniony. Ten katalog zawiera niezależne
weryfikatory i rejestr wyników przyjętych do repozytorium.

## Najważniejszy przyjęty wynik

Dla zadanych **dokładnych amplitud dziesiętnych** obie funkcje fazowe mają
**dokładnie 60 punktów krytycznych na całym torusie fazowym**:

- funkcja kwartowa K4;
- pełna funkcja log-mgf.

W obu przypadkach liczby punktów o indeksach ujemnych hesjanu 0, 1, 2, 3
wynoszą odpowiednio **12, 24, 18, 6**. To istotne wzmocnienie poprzedniego
wyniku „co najmniej 60 lokalnie certyfikowanych pierwiastków”: tym razem
sprawdzono również wykluczenie dodatkowych pierwiastków poza ich otoczeniami.

Nie jest to twierdzenie dla dowolnych amplitud, całego modelu X7 ani o fizycznym
wyborze stanu. Nie dowiedziono globalnej homotopii między K4 a pełną funkcją.

## Co sprawdzono niezależnie

| Warstwa | Zakres kontroli |
|---|---|
| Pochodzenie | Wszystkie 740 wpisów manifestu pakietu; trzy historyczne wejścia rozpoznane lokalnie po hashach. |
| Współczynniki | Nowa rekonstrukcja szeregu logarytmu inną rekurencją niż generator dostarczony w handoffie. |
| Błąd przybliżenia | Analityczne oszacowanie Cauchy’ego, dodatnia dziedzina analityczności i suma odrzuconych rezonansów K16/K20. |
| Lokalne pierwiastki | Ponowna certyfikacja 60+60 pierwiastków i indeksów hesjanu na podstawie poprawionego modelu z przedziałowymi stałymi. |
| Większe otoczenia | Nowe sprawdzenie kontrakcji/iniektywności w zadanych otoczeniach; nie tylko odczyt zapisanego `PASS`. |
| Wykluczenia K4 | Wszystkie 27 272 komórki wykluczające dodatkowe pierwiastki. |
| Wykluczenia pełnej funkcji | Wszystkie 54 341 komórek K16 oraz 25 292 komórki K20. |
| Pokrycie | Pełne podziały torusa i zawieranie 640/864 komórek w otoczeniach jednoznacznych pierwiastków. |
| Symetria | Dokładne działanie grupy zachowującej znak oraz przedziałowe przypisanie obrazów pierwiastków. |

Łącznie ponownie obliczono **106 905 nierówności wykluczających pierwiastki**.
Nie zastąpiono tej kontroli próbkowaniem 240 komórek ani wcześniejszymi flagami
ukończenia replay. Fazy są oceniane w znormalizowanych współrzędnych torusa;
testowana pochodna jest pochodną względem kąta φ, której zerowanie jest
równoważne zerowaniu pochodnej względem φ/(2π).

Przeszło również osiem nowych testów kontrolnych, m.in. odrzucenie brakującej
komórki pokrycia, niepełnego replay, otoczenia bez kontrakcji i kompresji
o niepełnej randze. Historycznych testów 119+7 nie liczono ponownie jako
nowych testów tego audytu.

## Wynik dla ograniczenia 4D

Przyjęty wynik jest **częściowy**, przy dokładnym progu `67/250`:

- 13 231 zaakceptowanych komórek;
- 5 432 komórki nierozstrzygnięte;
- około **63,6895053%** objętości zadanej zwartej dziedziny jest objęte dowodem;
- około **36,3104947%** pozostaje nierozstrzygnięte.

Sprawdzono geometrię podziałów dokładnym rachunkiem wymiernym oraz ponownie
obliczono nierówności na zaakceptowanych komórkach. Bazy używane w kompresji
macierzy wczytywano jako propozycje wymierne, po czym sprawdzano ich rangę
i dodatniość odpowiednich form. Nie polegano na zapisanych numerycznych
wartościach własnych lub samych flagach testów.

Powyższe procenty dotyczą **objętości współrzędnych zwartej dziedziny**,
nie prawdopodobieństwa fizycznego ani „procentu pewności twierdzenia”.
Globalne Target P i ostrzejsze Target S nadal pozostają otwarte.

## Uwagi metodologiczne i przenośność

1. `verify.py` oraz `portable_verify.py` nie są samodzielnym pełnym replay
   wszystkich dowodów. Weryfikator przenośny sprawdza m.in. zapisane flagi
   ukończenia i deterministyczne próbki. W tym audycie sprawdzono pełną warstwę
   matematyczną niezależnymi programami.
2. Część dostarczonych generatorów ma bezwzględne ścieżki `/mnt/data/...`.
   Nowe weryfikatory w tym katalogu korzystają ze ścieżek względnych repo.
3. Manifest obejmuje cztery pliki `.pytest_cache`. Pierwsza próba na kopii
   oczyszczonej z cache nie przeszła kontroli manifestu z tego właśnie powodu;
   nie oznacza to uszkodzenia oryginalnego pakietu. Wszystkie jego wpisy pasują.
4. Próba oryginalnego `portable_verify.py` na dokładnej kopii przekroczyła
   wewnętrzny limit 120 sekund podczas kontroli geometrii. Nie uznano tego
   za kontrprzykład ani błąd twierdzenia. Niezależne kontrole ukończono osobno.
5. Historyczne ścieżki trzech wejść są niedostępne pod oryginalnymi nazwami,
   ale ich zawartość rozpoznano w repozytorium i potwierdzono zgodność SHA-256.
6. Symetryczne przypisanie pierwiastków w dostarczonym skrypcie wykorzystywało
   numeryczne odległości. Tutaj dodatkowo sprawdzono przedziałowe zawieranie
   obrazu małego pudełka w otoczeniu z udowodnioną jednoznacznością.
7. Zachowano oddzielnie numeryczne propozycje, dowody przedziałowe i ograniczenia
   stosowalności. Nie rozszerzano wyniku fazowego na amplitudy odpowiadające
   nieznanemu dokładnemu punktowi koegzystencji.

## Przyjęte badania i artefakty

- [Dokładne twierdzenia i zakresy](ACCEPTED_RESULTS.md).
- [Niezależna rekonstrukcja K16/K20](surrogate_reconstruction.json).
- [Pierwiastki i otoczenia jednoznaczności](collar_replay.json).
- [Kompletność podziałów fazowych](phase_geometry.json).
- [Wykluczenia K4](leaf_replay_quartic.json), [K16](leaf_replay_16.json),
  [K20](leaf_replay_20.json).
- [Podgrupa symetrii i orbity](symmetry_replay.json).
- [Geometria częściowego pokrycia 4D](partial_geometry.json)
  i [replay zaakceptowanych komórek](partial_leaf_replay.json).
- [Pochodzenie i hashe](inventory.json).

Polecenia pełnego replay, z katalogu głównego repozytorium:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/review.py inventory
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/review.py reconstruct_surrogates
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/review.py collars
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/review.py phase_geometry
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/review.py leaves --kind quartic
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/review.py leaves --kind 16
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/review.py leaves --kind 20
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/review.py symmetry
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/partial_cover.py
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7n_review/finalize.py
```

[verification.json](verification.json) podsumowuje zakończone replay i hashe.
Ostatnie polecenie kontroluje ten stan i uruchamia osiem testów; samo nie
zastępuje wcześniejszych poleceń przeliczających dowody.

Skrypty zapisują wyniki audytu tylko w tym katalogu. Oryginalnych materiałów
nie poprawiono „w miejscu”, aby nie zniszczyć historii i manifestów.
Nie wykonano commita ani nowej kampanii mającej zamknąć pozostałe obszary 4D.
