# Audyt i integracja handoffu FIN — przed i po Discord

Data: 2026-09-13.

## Wynik audytu

Pakiet zawiera poprawne wyniki matematyczne i użyteczne wyniki numeryczne,
ale jego podsumowanie zawyża stopień domknięcia części badań. Przyjęto
wyniki strukturalne w jawnych zakresach, odtworzono główne obliczenia lokalne
i zbudowano dwa niezależne certyfikaty dla wyróżnionej powierzchni modelu.
Nie przyjęto twierdzenia o globalnym domknięciu modelu rank-7.

Co ważniejsze, znaleziono dokładny kontrprzykład dla przeniesienia postulatu
„wszędzie najwyżej jeden ujemny kierunek hesjanu” na wszystkie siedem
współrzędnych. Cztery amplitudy po ustaleniu faz i pełne siedem współrzędnych
to różne zadania. Kontrprzykład dotyczy hesjanu w dowolnym punkcie, nie
stanowi sam w sobie punktu stacjonarnego o indeksie dwa.

Źródła lokalne:

- [Oryginalny handoff](../FIN_research_artifacts_pre_and_post_Discord/FIN_full_chat_research_handoff_pre_and_post_Discord.md).
- [Manifest](../FIN_research_artifacts_pre_and_post_Discord/MANIFEST.txt):
  sprawdzono rozmiary i SHA-256 wszystkich 41 artefaktów: 28 CSV i 13 PNG.
- [Dowody i poprawione wzory](PROOF.md), [implementacja](research.py),
  [testy](test_research.py), [wyniki](results.json), [weryfikacja](verification.json).

Oryginalne tabele i wykresy pozostają niezmienione jako materiał źródłowy.
Ich obecność w repozytorium nie oznacza uznania wszystkich podpisów za dowody.
Wykresy nie były używane jako samodzielne świadectwa matematyczne.
Nie kontynuowano pobierania rozmowy internetowej po wskazaniu pakietu.

## Metodologia: co sprawdzono i czego brakuje

Przeczytano cały handoff (2577 wierszy), wszystkie tabele CSV i odpowiednie
aktualne ograniczenia w AGENTS.md, w tym ST293, ST344–ST361, ST389–ST410
oraz nowsze pakiety kwantowe. Porównano deklaracje z ich przesłankami,
zakresem przestrzeni stanów i faktycznie dołączonym materiałem.

Pakiet nie zawiera generatorów tabel, listy 60 pierwiastków, kompletnych
liści pokrycia Bernsteina, preconditionerów/pudeł Krawczyka ani pełnego
rachunku przedziałowego dla globalnej nierówności Isinga. Tabela z liczbą
nierozstrzygniętych pudeł nie jest certyfikatem ich rozstrzygnięcia.
Brak tych danych nie obala hipotez, lecz uniemożliwia przejęcie statusu
„computer-assisted theorem” bez odtworzenia dowodu.

Nowy kod niezależnie buduje operator z definicji, a nie z CSV. Dokładna
warstwa używa istniejących wymiernych przedziałów spektrum strict i
arytmetyki Fraction. Warstwa numeryczna używa NumPy/SciPy. Dodatnie
współczynniki Bernsteina są sprawdzane na przedziałach; pozorne zero
na końcu przedziału zastąpiono dokładnym wyłączeniem czynnika `1-q`.
Testy numeryczne nie są przedstawiane jako dowody globalności.

## Najważniejsze poprawki

1. **Parametr rezolwenty (§35).** Poprawny jest
   `r=sech(sqrt(lambda3/6)*s3)`, nie `exp(-sqrt(lambda3/6)*s3)`.
   Dla s3≈1.703819057 poprawny wzór i bezpośredni rachunek macierzowy dają
   S≈0.057549460989; błędne podstawienie daje około 0.180520227.
   Tabela minimum była zgodna z parametryzacją sech, nie z tekstem handoffu.
2. **Pełny model a cztery amplitudy (§20, §39, §47).** Dla
   `h_j=2 cos(pi*j/2)` pełny hesjan 7D ma co najmniej dwa ujemne kierunki
   dla każdego g≥3.7. Przy g≈3.718344898 są to numerycznie około -0.06125315
   i -0.06125315. Dowód wymierny znajduje się w PROOF.md §7.
3. **Brakujący dowód brzegowy (§27, §33, §39).** Nie jest poprawne uznanie,
   że jedyną pozostałą luką jest obszar poza wyróżnioną powierzchnią:
   globalny certyfikat Isinga i wynik zależny dla W_par również nie zostały
   dostarczone w odtwarzalnej postaci. Nowy dowód powierzchni nie zamyka tych luk.
4. **Warunki rezolwenty (§34, §37).** Wymagana jest odwracalność
   `sigma I-W_par`, a dla wskazanego testu znaku dokładnie jedna wartość
   własna W_par powyżej sigma. Ogólna poprawna redukcja używa bezwładności
   hesjanu i dodatniego mianownika Schura; przypadki równości są osobne.
5. **Szum walkerów (§2.3).** Macierz A/(6N) jest wartością przy rozkładzie
   jednostajnym, nie uniwersalną warunkową kowariancją w każdym stanie.
   Proces OU to liniowy opis fluktuacji; dokładny proces jest skokowy.
6. **Granice discord (§17, §44).** Wymierne dolne ograniczenia nie są równe
   dokładnym transcendentalnym lukom spektralnym. Zachowano dotychczasowy
   certyfikowany d0. Analogiczny współczynnik w postulowanej dynamice aktywnej
   nie dowodzi mechanizmu przyczynowego.
7. **Nazwa punktu odcięcia (§16, §46).** Opis 22 składników i ranku 132
   odpowiada `FIN_Separable_Stationarity_and_Necessary_Discord.tex`.
   Późniejszy raport `FIN_Discord_Robustness_and_Operational_Identifiability.tex`
   zawiera dalsze wyniki; nie wolno traktować streszczenia starszego raportu
   jako pełnego audytu późniejszego. Oba zachowano i uruchomiono ich regresje.

## Ledger obejmujący cały handoff

„Przyjęte” oznacza zakres opisany w PROOF.md, nie wszystkie zdania źródła.
„Numeryczne” nie oznacza certyfikacji lub globalnego wyczerpania.
„Archiwum” oznacza zachowanie dołączonych rezultatów bez promocji naukowej.

| Sekcje handoffu | Decyzja po audycie |
|---|---|
| 0–1 | Zachować rozdział przesłanek i statusów; spektrum odtwarzać z operatora i przedziałów, nie utożsamiać zaokrągleń z exact. |
| 2.1–2.2 | Przyjęte: rozkład 11/55, PSD, zerowa dywergencja, brak zależności od przecięć rysunku. Liczba 28 różnych dodatnich wartości własnych pozostaje numeryczna. |
| 2.3 | Przyjęte po korekcie: dokładny skokowy model walkerów, równowagowa wariancja 1/(12N), OU tylko jako opis liniowy. |
| 3 | Przyjęte: pasywność Schura; odtworzono rank 5, trace≈2.05511490123 i PSD operatora efektywnego. Minimalność realizacji/podział na 3 grupy biegunów nie jest tu nowym certyfikatem. |
| 3.1 | Model mean-field/Gibbsa pozostaje warunkową możliwością; brak pełnego mikroskopowego generatora w pakiecie. Nie źródło g. |
| 4 | Przyjęte i certyfikowane: rank 7 jako próg vertex/uniform przy g=4 w klasie D12, 0≤B≤A. Sygnatury rank 3 i generowanie częstotliwości sprawdzone skończenie. |
| 5 | Archiwum numerycznych kandydatów w drabinie ranków. „Pierwszy globalny competitor” nie jest dowiedziony; uwaga dotyczy także identyfikacji lokalnego crossing rank 11. |
| 6 | Przyjęte: dokładny dual 7D i nieliniowe domknięcie entropijne. Odtworzono ≈32.9075% mocy poza aktywną podprzestrzenią u kandydata. Liczby 4.720214 i 21.2% pozostają importowanymi wynikami optymalizacji ograniczonej, bez niezależnego certyfikatu. |
| 7 | Odtworzono lokalny crossing g≈3.71834489812038, pmax≈0.83636522665 i dodatni hesjan styczny. Fold≈3.5156447 pozostaje importowaną numeryką; spinodal 12/lambda6 jest tożsamością. Żadnego globalnego/przyrodniczego przejścia nie przyjęto. |
| 8 | Diagram 1→2→12 i promienie to wyniki poszukiwania kątowego; nie globalna klasyfikacja. Równanie niestabilności wymaga kontroli wszystkich kierunków; brak tej certyfikacji w dostarczonym pakiecie. |
| 9 | Rezonanse wynikają z dodawania częstotliwości mod 12; szczegółowy theorem 12 maksimów krajobrazu fazowego nie został tu certyfikowany. Nie przenosić minimum przy ustalonych amplitudach na pełną przestrzeń. |
| 10–11 | Zachować jako numeryczny katalog 60 znalezionych punktów i empiryczną zgodność obcięcia czwartego rzędu. Brak listy pierwiastków/pokrycia i jednostajnych granic C2; „exactly 60” nie jest twierdzeniem. |
| 12–14 | Archiwum: 95 znalezionych konkurentów, bariery radialne, bliskość promienia jednej cechy. Brak dowodu wyczerpania, braku maksimów lub dokładnej redukcji 1D. |
| 15 | Zachować aktualne ograniczenia repo; globalny theorem dla pełnego A przy g=4 nie transferuje do A7. |
| 16 | Istniejący wynik o separowalności i koniecznym discord, z poprawionym przypisaniem raportu. Nie liczyć ponownie jako nowego odkrycia. |
| 17 | Przyjęta algebraiczna analogia spektralna po rozdzieleniu dokładnych luk, dolnych enclosures i dostarczonego prawa wzrostu. Nie związek przyczynowy. |
| 18 | Przyjęta dokładna reprezentacja CRT; ogólna izotoniczność/covariance positivity nie została dowiedziona w tym audycie. Sama nazwa „ferromagnetyczny” nie płaci tego twierdzenia. |
| 19 | Odtworzono lokalny saddle i jeden ujemny kierunek w 4D. Globalna monotoniczna iteracja oraz wyczerpanie punktów stałych pozostają niecertyfikowane. |
| 20–21 | sigma≈0.267443244229 jest dokładnym kandydatem i certyfikowanym sufitem na opisanej powierzchni; nie globalnym sufitem pełnego 7D. Lokalny wyróżniony punkt brzegowy jest odtworzony. |
| 22 | Wskazane pojedyncze granice dużego pola opisują właściwe podpory; ogólne ścieżki ku nieskończoności wymagają osobnego pokrycia. Próbki i numeryczne szukanie gradientu nie są wyczerpaniem. |
| 23–24 | Zachować opisy nieudanych strategii jako wskazówki, nie gotowe kontrcertyfikaty. Wzór P'''=24sigma−6Tr(M) jest tożsamością; brak pełnego certyfikatu wskazanej kuli i implikacji znakowych. |
| 25–26 | Przyjęte: redukcja Isinga, degeneracja Y, trzy nierówności ilorazowe dla dodatnich p, dokładny determinant. Nie globalny curvature theorem. |
| 27 | NIE PRZYJĘTO jako computer-assisted theorem: 10 nierozstrzygniętych pudeł i opis lokalnego stożka nie zastępują replayowalnego zamknięcia. |
| 28 | sigma3≈0.200205788 pozostaje kandydatem globalnego sufitu, nie twierdzeniem. |
| 29 | Przyjęta dokładna dekompozycja całkowitej kowariancji. Trzeba zachować wspólne pola obu klas parzystości. |
| 30 | Przyjęto lokalne wzory pierwszego rzędu przy ustalonym tstar; ujemne znaki wynikają z certyfikowanych stałych. Nachylenie po globalnej reoptymalizacji pozostaje numeryczne. |
| 31–33 | Schemat Weyla jest warunkowo poprawny. Brak odtwarzalnego globalnego lematu Isinga i pokrycia dominującej masy: nie ogłaszać zamknięcia W_par. Tabele zachowane jako cele odtworzenia. |
| 34 | Przyjęto kryterium bezwładności z brakującymi przesłankami. Globalne minimum S≈0.05755 nie zostało dowiedzione. |
| 35 | POPRAWIONO sech; przyjęto nowy dokładny strict-interval certyfikat dodatniości na 1D powierzchni. Numeryczne minimum odtworzone, jego unikalność niecertyfikowana. |
| 36 | Archiwum numerycznych odległości od powierzchni i pochodnych; brak jednostajnej kontroli poza powierzchnią. |
| 37 | Przyjęta poprawiona ogólna równość liczb ujemnych kierunków po redukcji Schura, nie bezwarunkowa równoważność znaku S. |
| 38 | PRZYJĘTO: niezależny certyfikat dla całej powierzchni s4=s5=0, s3,s6≥0, już z przedziałami strict. Dokładnie opisana równość tylko na brzegu skompaktyfikowanym. |
| 39 | ODRZUCONO twierdzenie, że pozostała tylko jedna luka. Ponadto pełny everywhere-7D transfer jest obalony; restricted-4D i stationary-only to odrębne otwarte pytania. |
| 40 | Zachować historię nieudanych tras. Nowe twierdzenie o ich obaleniu wymaga zapisanych świadków; brakujące świadectwa nie są wytwarzane z opisów. |
| 41–42 | Dopuszczalna warunkowa interpretacja i ograniczenia ontologiczne. Bez źródła gainu, jednostek, selektora, laboratorium lub ToE. |
| 43 | Priorytety skorygowane: część podniesienia przedziałowego wykonana; nadal brakuje brzegowego cover, kontroli poza powierzchnią i kontroli faz. Nie kontynuować automatycznie nowej kampanii globalnej. |
| 44 | Stałe to mieszanka tożsamości, przedziałów i kandydatów numerycznych; statusy należy brać z tego audytu. |
| 45–46 | Kompletność 41 plików potwierdzona hashami. Dodatkowe pierwotne wklejone teksty i historyczny AGENTS(1).md nie należą do dostarczonych artefaktów; nie udajemy audytu ich pełnej treści. |
| 47 | Podsumowanie zastąpione niniejszym werdyktem: istotne wyniki lokalne są realne, ale „near-closure” wymaga więcej niż samej kontroli s4/s5. |

## Co faktycznie włączono

- Odtwarzalną implementację operatora, obu dualnych przestrzeni współrzędnych,
  podziału parzystości, lokalnego crossing i saddle oraz struktur 11/55 i Schura.
- Dokładne certyfikaty przedziałowe progu budżetu oraz dwóch wyników
  powierzchniowych Bernsteina; zakresy i dowody w PROOF.md.
- Dokładny kontrprzykład blokujący nadinterpretację pełnego modelu 7D.
- Ledger wszystkich sekcji źródła, poprawki i jawne statusy otwartych hipotez.
- Nową sekcję AGENTS.md odsyłającą do tego audytu.

Weryfikacja:

```bash
python3 fin_handoff_audit/verify.py --record
python3 fin_handoff_audit/verify.py
```

Pierwsze polecenie zapisuje nowe wygenerowane ledgery, drugie sprawdza ich
dokładną warstwę i uruchamia regresje. Zestaw obejmuje 19 nowych testów
oraz 56 testów trzech wcześniejszych pakietów kwantowych. Te regresje
sprawdzają zgodność istniejących implementacji; nie stanowią niezależnego
odtworzenia każdej historycznej argumentacji tych pakietów. Importowane
globalne pokrycia i katalog 60 pierwiastków nie są odtwarzane, ponieważ
ich nie dostarczono. Nie wykonano commitów, nie wygenerowano PDF,
nie wznowiono zatrzymanych ciężkich kampanii obliczeniowych.
