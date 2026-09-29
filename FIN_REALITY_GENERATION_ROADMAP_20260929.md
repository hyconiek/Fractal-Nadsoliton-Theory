# FIN — roadmap sprawdzenia hipotezy generowania struktur fizycznych

Data: 2026-09-29.
Stan odniesienia: cf9b6607.
Zakres nowej analizy: fin_physical_bridge_campaign_01 i fin_physical_bridge_campaign_03; wcześniejsze guardraile 327–337 pozostają w mocy.

To aktualizacja kierunku po FIN_PHYSICS_BRIDGE_MASTER_ROADMAP_20260929.md.
Nie ponawia wykonanych PHYS-001–007. Nie jest deklaracją, że FIN generuje
Wszechświat. Nie zmienia AGENTS.md, protokołu 333, statusu 335/337 ani
nie uruchamia eksperymentu sprzętowego.

## 1. DIAGNOSIS — co zmieniły najnowsze pliki

### Werdykt

FIN ma obecnie:
- dokładny model-class bridge do generalized matrix Curie–Weiss;
- zależne od A7 statyczne fingerprinty;
- inżynierski projekt testu dwóch q=12 węzłów;
- nadal NO_NEW_SOURCE dla konkretnej A7.

Nie ma jeszcze dowodu, że FIN wytwarza z własnych zasad:
A7, fizyczną kompozycję podukładów, przestrzeń, fundamentalną kinetykę,
QM, grawitację ani rzeczywistą ontologię świata.

Najbliższym uzasadnionym celem generatywnym jest:
ustalić, jakie nowe prawa efektywne powstają z obecnego modelu po
eliminacji części stopni swobody, jak kontrolować ich błędy i co
odróżnia ten mechanizm od standardowej klasy mean-field.

### Aktualne wyniki

| Pakiet | Co wynika | Czego nie wynika |
|---|---|---|
| PHYS-001–004 | Dokładne mapowanie i finite-N static fingerprint; oddzielenie równowagi od kinetyki | fizyczna konieczność A7 |
| PHYS-005 | Ledger NO_NEW_SOURCE; Gaussian factorization nie wybiera konkretnego kernelu | prawo źródłowe z samej PSD/symetrii |
| PHYS-006–007 | Projekt realizacji i zamrożony kontrast N=2 | laboratoryjna walidacja |
| PHYS-011 | Dwa q=12 węzły mają dokładnie ten sam idealny histogram różnicy co wcześniejszy test conditional | osiągnięcie deklarowanej dokładności przez fizyczny hardware |
| PHYS-012 | Nie można dziedziczyć fidelity 0.003 z ogólnych wyników istniejącej maszyny Pottsa | pełny rozdział biasu sprzętu od błędu próbkowania skorelowanych danych |
| PHYS-013 | Nominalny 144-state CTMC kompilacji daje TV około 0.0004127833 | zmierzony TV rzeczywistej konfiguracji FIN |
| PHYS-014 | READY_FOR_REAL_Q12_PRECALIBRATION_ONLY | otwarcie validation, zakup lub kontakt z laboratorium |

Sprawdzono manifesty obu nowych kampanii. Odtworzono mały PHYS-013
po zmianie wyłącznie ścieżki importu w strumieniu odczytu:
TV=0.0004127833002122317. To replay modelu z zapisanymi parametrami
kalibracji, nie niezależne odtworzenie fizycznych pomiarów.

Plik ctmc_replay.py wymaga poprawienia przenośności importu:
odwołuje się do nieobecnego root fin_core.py lub historycznego /mnt/data.
Właściwy moduł i jego hash muszą być związane z finalnym replay.

### Ważne uzupełnienia metrologiczne

1. W PHYS-012 liczba zapisanych próbek nie jest automatycznie liczbą
   niezależnych próbek. Oszacowanie 0.5 sqrt(K/n) wymaga odpowiednich
   założeń o próbkowaniu. Bez analizy autokorelacji/ESS nie można samym
   tym porównaniem wykluczyć wyjaśnienia części rozbieżności przez sampling.
   Wniosek „nie dziedziczyć fidelity” pozostaje prawidłowy.
2. Nominalny CTMC musi być uzupełniony o uncertainty temperatur kanałów,
   dryft, cross-talk, opóźnienia, invalid states i serial correlation.
   Błąd 0.000413 nie jest górnym ograniczeniem tych niezmierzonych efektów.
3. Zgodność stacjonarnego histogramu nie dowodzi realizacji pełnej
   leave-one-out dynamiki. Statics i dynamics mają osobne bramki.
4. Referencyjna publikacja potwierdza istnienie SPAD Ising/Potts hardware
   i podaje źródła 16/144-circuit data. Nie opisuje walidacji FIN q=12.
   Źródła:
   https://www.nature.com/articles/s41928-023-01065-0
   https://github.com/ucsb-biomimetic/ising-potts-with-SPADs-discrete
   Konkretne raw calibration files należy związać z URL, revision i hash,
   a nie tylko z nazwą i nieadresowalnym blob SHA.

## 2. Co znaczy „FIN generuje rzeczywistość”

Należy rozdzielić cztery poziomy.

**G0 — generowanie trajektorii.**
Program symuluje zadany proces. To już możliwe. Nie dowodzi fizyki.

**G1 — generowanie obserwowalnej fizyki efektywnej.**
Jedno jawne prawo mikroskopowe przewiduje fazy, korelacje, odpowiedzi
i procesy poza kalibracją. FIN ma tu matematyczne wyniki i kandydat
realizacji, ale brak ukończonego nowego testu fizycznego.

**G2 — generowanie struktury teorii.**
Nie wpisuje się docelowej A7, konkretnego grafu ani kolejnego generatora
osobno na każdym poziomie. Wynikają one z mniejszego, wcześniej
określonego zestawu zasad i mają niezależne testowalne konsekwencje.
Tego obecny FIN jeszcze nie wykazał.

**G3 — wyprowadzenie znanej rzeczywistości fundamentalnej.**
Lokalność, czasoprzestrzeń, pola/materia i ewentualne QM/GR wynikają
z wcześniejszej struktury i przechodzą własne testy.
To daleki cel warunkowy, nie następny gotowy etap.

Celem następnej kampanii jest sprawdzić przejście G1→fragment G2.
Nie wolno raportować sukcesu G0 lub realizacji zaprogramowanego A7 jako G3.

Zasada przeciw kołowości:
obiekt nie jest „wygenerowany”, jeśli został dostarczony równoważnie
w macierzy, funkcji kosztu, wyborze generatorów, historii przyszłości,
warunkach początkowych albo aparacie.

Przyjęta ontologia projektu pozostaje bez zmian: nadsoliton jest
pierwotną informacją w stanie solitonicznym. Nie dodaje się warstwy
informacyjnej pod nim. Nie jest to jednak wynik empiryczny tej kampanii.

## 3. CLOSEST PHYSICS i nadal obowiązujące mapowanie

Obecny rdzeń:

H = -J_c/(2N) Σ_(a,b) A7[sigma_a,sigma_b] - h n0,
g=beta_th J_c, theta=beta_th h.

To generalized q=12 Curie–Weiss / discrete vector-spin mean-field.
W urządzeniu programowanym H może być jedynie Hamiltonianem celu,
nie jego rzeczywistą energią elektryczną.

| FIN | Rola kontrolowana | Nie wolno automatycznie wnioskować |
|---|---|---|
| 12 etykiet | 12 stanów wewnętrznych | 12 cząstek lub punktów przestrzeni |
| A7 | macierz oddziaływań modelu | konieczny kernel natury |
| g | beta J_c albo kalibrowany gain | czas kosmologiczny |
| theta/kappa | zadany bias przygotowania | fundamentalna stała |
| N | liczba kopii w modelu/realizacji | liczba bitów Wszechświata |
| rho | wolne tempo modelu razy kalibrowana częstotliwość | samodzielnie wyprowadzona sekunda |
| stan FPGA/SPAD | realizacja konkretnego procesu | wyprowadzona ontologia świata |
| hidden-variable elimination | nowe oddziaływania i pamięć efektywna | nowe fundamentalne siły bez dodatkowego mostu |

Najważniejszy aktualny diagram sprzętowy jest:
A7 jako wejście → skompilowane wagi → proces → histogram.
Nie jest to diagram:
fizyka urządzenia → niezależnie wyprowadzona A7.

## 4. Nowy kierunek generatywny: domknięcie po zmianie opisu

### Dlaczego to nowy problem, a nie kolejny census siodeł

Jeżeli FIN ma prowadzić do teorii wielu skal, trzeba wiedzieć, jakie
oddziaływania powstają po usunięciu części kopii. Nie ma powodu,
aby rodzina „jedna A7 i jeden gain” pozostawała dokładnie zamknięta.

Dla N+1 oznakowanych kopii z zerowym polem i stałą przekątną A7,
po wyeliminowaniu ostatniej kopii względna log-waga pozostałych N wynosi,
z dokładnością do stałej:

log W_eff(n)
= g/[2(N+1)] n^T A7 n
  + log Σ_j exp{g/(N+1) (A7 n)_j}.

To dokładna tożsamość dla wagi konfiguracji oznakowanych.
Jeżeli używa się liczebności, osobno trzeba uwzględnić znany czynnik
multinomialny; nie traktować go jako nowego oddziaływania.

Drugi składnik generuje ogólnie operatory wyższych rzędów.
Nie należy go zastępować samym dopasowaniem nowego g.

### Mały świadek kierunku

Dla przejścia cztery kopie → trzy kopie i już używanego gainu
g=5.145228719489142, sprawdzono osiem konfiguracji:

sigma1 ∈ {0,1}, sigma2 ∈ {0,2}, sigma3 ∈ {0,4}.

Naprzemienna różnica trzeciego rzędu ich log-wag wyniosła numerycznie:

Delta_123 log W_eff ≈ -0.27127058458855036.

Każda log-waga składająca się wyłącznie ze stałej, pól jednociałowych
i oddziaływań parowych ma taką różnicę równą zero.

To silny NUMERYCZNY WITNESS kandydatury interakcji trójciałowej;
należy spisać dowód i zapłacić arithmetic/input error przed
oznaczeniem go jako formalny certyfikat. Nie otwarto nowego gainu
ani dużej kampanii.

Wniosek:
po redukcji FIN trzeba śledzić nowe operatory, a nie kopiować A7
pod inną nazwą. Nie jest to refutacja FIN: to normalna możliwość
w teorii efektywnej. Nie jest też samo w sobie FIN-specific —
kontrmodele muszą przejść tę samą procedurę.

### Co należy sprawdzić dalej

- które operatory są rzeczywiście generowane;
- które są potrzebne do obserwowalnych przewidywań;
- czy ich współczynniki wynikają z jednego mikromodelu bez nowych fitów;
- czy kolejne redukcje są zgodne przy jawnej kontroli błędu;
- czy FIN różni się tu od zwykłych kontrmodeli tej samej klasy.

To jest marginalizacja podukładów mean-field, nie przestrzenny Wilson RG.
Sama operacja nie tworzy lokalności ani fizycznych długości.

## 5. CRITICAL MISSING BRIDGES

### Źródło oddziaływania

Istniejący NO_NEW_SOURCE pozostaje w mocy. D12, PSD, trace i rank
nie wybierają konkretnej A7.

Najbardziej konkretny niezależny source-test pozostaje taki:
z rzeczywistej, wcześniej określonej fizyki mediatorów otrzymać K i c_j,
a następnie wyprowadzić A_eff ∝ c_i^T K^-1 c_j.

Nie wolno:
- wybrać c_j przez faktoryzację gotowej A7 i nazwać tego źródłem;
- dobrać liczby modów i ich stiffness po obejrzeniu FIN fingerprint;
- uznać fizycznego pomiaru zaprogramowanych współczynników za dowód,
  że natura niezależnie wybiera ten kernel.

Sukces takiego testu najpierw oznaczałby „FIN jako opis efektywny
konkretnego mechanizmu”. Jeśli parent model zakłada przestrzeń i
czas, nie wyprowadzono ich z FIN.

### Autonomia i warunki początkowe

Aktualne przygotowanie korzysta z deep seed, biasu i harmonogramu.
To jawne wejścia. Aby mówić o autonomicznym generatorze, trzeba
określić stan sterownika, środowiska i źródła zdarzeń.

Natural extension dowodzi możliwości odwracalnego opisu historii,
nie fizycznego prawa, które wybiera tę historię.
Nie powtarzać prób absolutnej kompletności z tych samych projected data
po no-go 295.

### Kompozycja i lokalność

Trzeba wyprowadzić:
jednoczesne operacyjne podukłady → prawo ich wspólnej dynamiki →
kontrolowane rozchodzenie się wpływu.

All-to-all mean-field nie staje się lokalną przestrzenią po
przerysowaniu grafu. Geometria stanów i geometria jednoczesnych
podukładów to różne obiekty.

### Skala, pola i materia

Po źródle kompozycji potrzebne są:
- spójna rodzina zwiększanych układów;
- określone transformacje skali;
- stabilne obserwable i kontrola pamięci;
- wyprowadzone prawa zachowania;
- efektywne pola jako zachowane wolne zmienne;
- dopiero potem test stabilnych propagujących wzorców.

Gap relaksacyjny nie jest automatycznie masą cząstki.
Basin multiplicity nie jest liczbą rodzin materii.

### Późne bramki QM/GR

QM wymaga osobnego operacyjnego mostu do interferencji, kompozycji,
niekomutujących pomiarów i reguły prawdopodobieństw.
Zespolony Fourier lub unitary representation permutacji nie wystarcza.

GR wymaga co najmniej dynamiki efektywnej geometrii, backreaction,
praw zachowania i właściwej struktury kauzalnej.
Stanowo zależna macierz przewodności nie jest jeszcze grawitacją.

Nie zlecać tych etapów przed wcześniejszymi bramkami.

## 6. MASTER ROADMAP

~~~text
AKTUALNY MODEL FIN
  |
  +--> PRZYJĘTY MODEL-CLASS BRIDGE
  |       -> dwa q12 węzły / rzeczywista prekalibracja
  |       -> niezależny test wobec kontrmodelu
  |       -> walidacja REALIZACJI, nie źródła A7
  |
  +--> GENERATYWNOŚĆ OPISU EFEKTYWNEGO
          -> dokładne marginalizacje
          -> wygenerowane operatory + kontrola truncation
          -> porównanie z kontrmodelami / uniwersalność
          -> CHECKPOINT 1
                 |
                 +-- brak source / brak swoistej predykcji
                 |      -> STOP fundamental-source inference
                 |
                 +-- nowy niezależny source object
                        -> wyprowadzenie kernelu/rodziny bez refit
                        -> autonomia + kompozycja podukładów
                        -> operacyjna lokalność
                        -> scaling / efektywne pola
                        -> przewidywanie spoza kalibracji
                        -> CHECKPOINT 2
                        -> dopiero warunkowo quantum/gravity bridges
~~~

Wybrany najbliższy kierunek:
kontrolowana generacja praw efektywnych przez marginalizację,
z równoległą, NIEZAMKNIĘTĄ bramką niezależnego źródła A7.

Nie obiecuje to generowania Wszechświata. Umożliwia tani test
koniecznego elementu takiej ambicji: czy kolejne poziomy opisu wynikają
z jednego modelu, czy są dopisywane ręcznie.

## 7. Dokładne wejścia nowej kampanii

| Alias | Ścieżka |
|---|---|
| A | AGENTS.md |
| B | FIN_PHYSICS_BRIDGE_MASTER_ROADMAP_20260929.md |
| C | fin_physical_bridge_campaign_01/HANDOFF.md |
| D | fin_physical_bridge_campaign_01/PHYS-001/DICTIONARY.md |
| E | fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py |
| F | fin_physical_bridge_campaign_01/PHYS-003/PROOF.md |
| G | fin_physical_bridge_campaign_01/PHYS-004/COUNTERMODELS.json |
| H | fin_physical_bridge_campaign_01/PHYS-004/CALIBRATION_CONTRACT.json |
| I | fin_physical_bridge_campaign_01/PHYS-005/SOURCE_GRAPH.md |
| J | fin_physical_bridge_campaign_01/PHYS-005/SOURCE_PREMISES.json |
| K | fin_physical_bridge_campaign_03/HANDOFF.md |
| L | fin_physical_bridge_campaign_03/PHYS-013/ctmc_replay.py |
| M | fin_physical_bridge_campaign_03/PHYS-013/RESULTS.json |
| N | fin_physical_bridge_campaign_03/PHYS-014/REVISED_PREREGISTRATION.md |
| O | fin_research_327_current_review/INTAKE_20260928.md |
| P | FIN_RESEARCH_AFTER_256_294_FULL_HANDOFF_295_299_20260927/01_REPORTS/COMPLETE_SYSTEM_CRITERION_295.md |

Nie potrzebować ponownie całej historycznej kampanii.
Brakujący input oznaczyć OPEN_INPUT, nie odtwarzać go z nagłówka raportu.
Nowa zawartość archiwum/rozpakowania FIN son nie stanowi tu przyjętego
wejścia: audyt jego archiwalnej treści nie został wykonany. Wcześniejsza
ocena dwóch jawnych skryptów pozostaje w dokumencie B.

## 8. NEXT CAMPAIGN — zadania wykonawcze

To plan, nie wykonane badania. Każdy task eksportuje REPORT.md,
RESULTS.json, NONCONCLUSIONS.md, REPLAY.md i MANIFEST.sha256.
Rozdzielać wykonanie od naukowego PASS. Nie poprawiać AGENTS.md automatycznie.

### GEN-001 — EXACT-MARGINAL-GENERATION

**Klasa / priorytet:** [B], P0.

**Cel:** jakie interakcje powstają po usunięciu jednej kopii?

**Dlaczego:** testuje strzałkę mikromodel → następny opis efektywny.

**Wejścia:** A, D, E; dokładna formuła z sekcji 4.

**Dependencies:** istniejący PHYS-001; nie ponawiać jego dowodu.

**Kroki:**
1. Wyprowadź marginalizację oznakowanych kopii dla N+1→N.
2. Oddziel multinomial degeneracy przy przejściu do count states.
3. Odtwórz świadek ośmiu konfiguracji dla 4→3 przy istniejącym g.
4. Spisz zerowanie trzeciej różnicy dla dowolnego pairwise Hamiltonianu.
5. Certyfikuj znak świadka z enclosures wejść i exp/log.
6. Sprawdź, czy 3→2 zachowuje Fourier support pierwotnej A7.

**Wynik:** dokładna tożsamość, minimalny świadek i opis operatorów, które trzeba zachować.

**PASS:** wyprowadzony effective log-weight i rozstrzygnięty pairwise-closure test.

**Kill-test:** jeżeli family A7+g nie jest zamknięta, zakazać jej kopiowania między poziomami bez error bound. Nie oznacza to automatycznej refutacji modelu.

**Forbidden shortcuts:** nowy g fit jako zastępstwo generated terms; pomylenie log-miary z microscopic force; global Γ=B4 jako nieopłacona przesłanka.

**Po PASS:** GEN-002.
**Po FAIL wykonania:** popraw mały checker; jeśli potrzeba nowej definicji modelu, wróć do [S].

### GEN-002 — MINIMAL-DERIVED-EFFECTIVE-ACTION

**Klasa / priorytet:** [B], P0.

**Cel:** czy wygenerowane oddziaływania można zachować w kontrolowanym, niewielkim opisie?

**Dlaczego:** następna skala nie może być nową niezależną tabelą fitowaną do wyniku.

**Wejścia:** GEN-001, D, E, F.

**Dependencies:** GEN-001 zakończone z opłaconym zakresem.

**Kroki:**
1. Użyj jednoznacznego gauge: stała + pola + pair terms + residual higher-body, np. rozkład względem uniform reference.
2. Zachowaj permutation symmetry kopii i D12 covariance etykiet.
3. Wyprowadź współczynniki z marginalizacji, nie z dopasowania obserwacji.
4. Zbuduj kolejno pair-only oraz najmniejszy model z wygenerowanym członem wyższym.
5. Ogranicz błąd log-wagi i przełóż go na błąd znormalizowanego rozkładu.
6. Porównaj redukcję dwuetapową z bezpośrednią marginalizacją; opłać truncation.

**Wynik:** operator ledger i jedna kontrolowana reguła redukcji.

**PASS:** error bound oraz brak nowego współczynnika dopasowanego oddzielnie po każdym kroku.

**Kill-test:** liczba koniecznych swobodnych fitów rośnie bez kontroli albo truncation nie ma użytecznego błędu.

**Forbidden shortcuts:** nazywanie marginalizacji przestrzennym RG; zachowanie starego ranku przez nieopłaconą projekcję.

**Po PASS:** GEN-003.
**Po FAIL:** zachowaj dokładną teorię z większą rodziną operatorów i skieruj wybór architektury do [S]; nie maskuj błędu.

### GEN-003 — GENERATIVE-SPECIFICITY-VS-UNIVERSALITY

**Klasa / priorytet:** [A], P0.

**Cel:** czy zachowanie przy redukcji wyróżnia FIN czy całą standardową klasę?

**Dlaczego:** wygenerowane interakcje wielociałowe same nie są unikatową fizyką FIN.

**Wejścia:** GEN-001/002, G, H, F.

**Dependencies:** GEN-001/002 PASS.

**Kroki:**
1. Użyj wcześniej zamrożonych full-Potts, flat-P7 i perturbation controls.
2. Zastosuj dokładnie tę samą mapę marginalizacji i ten sam error budget.
3. Porównaj powstające operatory, symmetrie, fingerprinty i wrażliwość.
4. Oddziel model-specific coefficients od strukturalnie stabilnych rezultatów.
5. Nie dopasowuj ponownie skal dla każdego obserwowanego efektu.
6. Wyeksportuj tabelę: swoiste / wspólne / nierozstrzygnięte.

**Wynik:** GENERATIVITY_CONTRAST.json i propozycja jednego dalszego falsyfikowalnego pytania.

**PASS:** kompletny scoped comparison, również jeśli rezultat jest negatywny.

**Kill-test:** jeśli wszystkie istotne efekty są wspólne klasie, nie używać ich jako dowodu wyróżnienia FIN.

**Forbidden shortcuts:** asymptotic fixed point z dwóch małych redukcji; claim o fizycznej przestrzeni z label graph.

**Po PASS lub naukowym NO_GO:** CHECKPOINT 1, czyli GEN-008. Nie zaczynać GEN-004–007 automatycznie.

### GEN-004 — INDEPENDENT-SOURCE-INPUT-GATE

**Klasa / priorytet:** [B], P0 dla silniejszej interpretacji.

**Cel:** czy pojawił się nowy, rzeczywiście niezależny obiekt źródłowy dla A7?

**Dlaczego:** aktualny NO_NEW_SOURCE nie znika po zbudowaniu samplera.

**Wejścia:** I, J, wynik checkpointu 1 oraz NOWY wskazany przez architekta materiał fizyczny. Bez tego ostatniego nie uruchamiać.

**Dependencies:** osobna decyzja [S] po GEN-003; nowy input poza istniejącymi target couplings.

**Kroki:**
1. Zapisz jego pochodzenie, role fizyczne, skalę i kontrolowane parametry.
2. Jeśli jest to mediator, wyznacz K i c_j z niezależnych reguł/pomiarów.
3. Wyprowadź A_eff i zidentyfikuj, które liczby są nadal ręcznie zadane.
4. Zamroź selection rules, normalizację i tolerancje przed porównaniem z FIN.
5. Sprawdź cutoff, rank i trzy proporcje, nie tylko możliwość PSD factorization.
6. Przygotuj przewidywanie dla zmiany fizycznego parametru spoza kalibracji.

**Wynik:** SOURCE_CERTIFICATE albo NO_NEW_SOURCE/OPEN_INPUT.

**PASS:** kernel lub jego kontrolowana klasa wynika z niezależnego inputu i ma nowe przewidywanie.

**Kill-test:** c_j pochodzą z faktoryzacji gotowej A7; indywidualne eigenvalues są fitted do targetu; source wykorzystuje ten sam fingerprint co test.

**Forbidden shortcuts:** brakujące dane zastąpione wymyślonym pomiarem; nowa zasada fundamentalna wymyślona przez wykonawcę; restart wyczerpanych symmetry-only no-go.

**Po PASS:** GEN-005 i source-test część GEN-007.
**Po FAIL/OPEN_INPUT:** stop source lane. Możliwy jest nadal oddzielny engineered bridge, bez promocji fundamentalnej.

### GEN-005 — AUTONOMOUS-COMPOSITION-AND-LOCALITY-GATE

**Klasa / priorytet:** [B], P1; warunkowo.

**Cel:** co nowy parent law naprawdę generuje po złożeniu dwóch i trzech jednostek?

**Dlaczego:** obraz świata wymaga kompozycji, nie tylko alternatywnych stanów jednej komórki.

**Wejścia:** GEN-004 PASS, zaakceptowana parent dynamics, P, aktualne guardraile.

**Dependencies:** niezależny source object albo jawnie zatwierdzona przez [S] nowa hipoteza modelowa; oznaczyć, który przypadek zachodzi.

**Kroki:**
1. Zapisz pełny state/controller/environment ledger i dozwolone operacje.
2. Wyprowadź joint law dwóch, następnie trzech jednostek z tego samego parentu.
3. Nie dodawaj nowej parowej macierzy, zegara ani grafu osobno przy każdej kompozycji.
4. Wyznacz wpływ interwencji, zachowane wielkości i pamięć po redukcji.
5. Sprawdź, czy lokalność jest wynikiem, czy założeniem parentu.
6. Jeżeli parent używa już przestrzeni/czasu, zaznacz te obiekty jako INPUT, nie GENERATED.

**Wynik:** composition theorem albo ograniczony kontrprzykład oraz tabela inputs/outputs.

**PASS:** operacyjnie spójna kompozycja z kontrolowaną pamięcią i bez ukrytego refit.

**Kill-test:** aktualny all-to-all model nadal daje pełny influence graph; nazwanie go przestrzenią nie przechodzi bramki.

**Forbidden shortcuts:** history slot=spatial site; conditioning on controller usuwa jego koszt; reversible gate=wyprowadzony autonomiczny zegar.

**Po PASS:** GEN-007, a plan scaling/fields dopiero po [S].
**Po FAIL:** zachowaj działający finite collective model; zamknij interpretację emergentnej lokalności.

### GEN-006 — ENGINEERED-BRIDGE-UNCERTAINTY-CLOSURE

**Klasa / priorytet:** [A], P0 w oddzielnej gałęzi pomiarowej.

**Cel:** czy nominalny projekt PHYS-013 może spełnić gate PHYS-014 przy rzeczywistych niepewnościach?

**Dlaczego:** simulated accuracy nie jest calibration certificate.

**Wejścia:** K–N, public source URL i zamrożony revision/data hashes, osobna zgoda człowieka przed realną akwizycją.

**Dependencies:** nie zależy od GEN-004; pozostaje testem inżynierskim.

**Kroki:**
1. Napraw portable import replay bez zmiany modelu.
2. Odtwórz pochodzenie 24 kanałów z pełnej tabeli, nie tylko wklejonych wierszy.
3. Przelicz perturbacje kalibracyjne i ich wpływ na stacjonarny histogram oraz currents.
4. Zdefiniuj uncertainty intervals, drift rules, invalid outcomes i próbki efektywne.
5. Nie używaj 4.9995 mln rekordów jako liczby iid bez uzasadnienia.
6. Przy braku realnych q12 records zakończ na READY_FOR_PRECALIBRATION.
7. Pomiar sprzętowy rozpocznij tylko po odrębnej autoryzacji, bez retuning validation.

**Wynik:** audytowalny calibration/error contract; później ewentualnie niezależny rekord fizyczny.

**PASS:** wszystkie składniki 0.003 są opłacone dla deklarowanego observable, albo jawnie zachowany status DESIGN_ONLY.

**Kill-test:** robust uncertainty region przecina kontrmodel; nominalne 0.000413 nie usuwa tego problemu.

**Forbidden shortcuts:** source claim z programmed weights; odrzucanie invalid events; Bayes conditional design jako laboratoryjna pewność.

**Po PASS:** empirical checkpoint [S], bez promocji źródła A7.
**Po FAIL:** popraw tylko niezależnie kalibrowalny element lub zatrzymaj platformę dla tej dokładności.

### GEN-007 — GENERATED-PREDICTION-OUTSIDE-CALIBRATION

**Klasa / priorytet:** [A]/[B], P1.

**Cel:** czy wyprowadzone prawo przewiduje nowy wynik fizyczny, a nie tylko odtwarza dane wejściowe?

**Dlaczego:** to właściwa bramka od mechanizmu generatywnego do fizyki.

**Wejścia:** GEN-004 i ewentualnie GEN-005, jawny physical input/output map,
wyniki calibration lane; bez nich tylko projekt i OPEN_EXTERNAL.

**Dependencies:** source claim wymaga GEN-004 PASS. Sam GEN-006 licencjonuje tylko engineered validation.

**Kroki:**
1. Wybierz jedną nową zmianę kontrolowanego parametru fizycznego.
2. Oblicz kernel/response law z parentu, bez wpisywania nowej A7 po zmianie.
3. Zamroź predykcję, nuisance intervals i kontrmodele.
4. Oddziel scalar clock calibration od bezwymiarowego fingerprint.
5. Porównaj nowe dane dopiero po freeze.
6. Zapisz wpływ błędu source, redukcji, przygotowania i odczytu.

**Wynik:** pre-registered prediction i późniejszy wynik niezależnego testu.

**PASS:** nowa obserwacja mieści się w uprzednio wyznaczonej klasie, przy rozróżnialnej alternatywie i bez refit.

**Kill-test:** nowy punkt wymaga swobodnej zmiany kernelu, grafu, jednostek lub definicji obserwabli.

**Forbidden shortcuts:** własna symulacja jako dane przyrody; source i validation z tego samego fitu; wynik jednego modelu jako dowód generowania całego świata.

**Po PASS/FAIL:** GEN-008.
**Po braku danych:** OPEN_EXTERNAL i stop; nie zastępować danych narracją.

### GEN-008 — ARCHITECT GENERATIVITY CHECKPOINT

**Klasa / priorytet:** [S], P0.

**Cel:** określić najwyższy rzeczywiście osiągnięty poziom G0–G3.

**Wejścia:** rezultaty poprzednich tasków, obecny NO_NEW_SOURCE, claim DAG i nonconclusions.

**Dependencies:** CHECKPOINT 1 po GEN-001–003; kolejny po ewentualnych source/experimental gates.

**Kroki:**
1. Oznacz każdy obiekt: INPUT, DERIVED, FITTED, EMPIRICALLY_TESTED.
2. Oddziel standardowe coarse-graining effects od FIN-specific prediction.
3. Rozstrzygnij, czy nowy source jest niezależny czy tylko inną reprezentacją targetu.
4. Sprawdź, które cechy świata parent zakładał z góry.
5. Wybierz dokładnie jedną dalszą gałąź albo STOP.

**Wynik:** CONDITIONAL_EFFECTIVE_MODEL, ENGINEERED_PHYSICAL_REALIZATION,
INDEPENDENT_SOURCE_SUPPORTED, GENERATIVE_CANDIDATE albo STOP_SOURCE_LANE.

**PASS:** verdict z opłaconym zakresem i konkretnym następnym falsyfikowalnym pytaniem.
Żaden z tych statusów nie oznacza „FIN wygenerował Wszechświat”.

**Kill-test:** wszystkie istotne elementy są wejściami albo uniwersalnymi efektami,
a kolejny krok ma tylko ukryć brak source w nowym fitcie.

**Forbidden shortcuts:** uruchamianie QM/GR/particle-number campaigns bez wcześniejszych gates.

**Po PASS:** osobny ograniczony plan następnego etapu.
**Po FAIL:** zachować użyteczną teorię efektywną, zatrzymać nieuzasadnioną interpretację fundamentalną.

## 9. CHECKPOINTY i instrukcja handoff

Pierwszy rzeczywisty powrót do architekta następuje po GEN-001–003.
Nie zaczynać source/composition/geometry campaigns automatycznie.
GEN-006 jest oddzielną gałęzią inżynierską; nie zastępuje source gate.

Gotowy tekst:

~~~text
Przeczytaj FIN_REALITY_GENERATION_ROADMAP_20260929.md i AGENTS.md.
Wykonaj wyłącznie GEN-001, GEN-002 i GEN-003, z istniejących danych
i przy już używanym g. To mała kampania exact marginalization;
maksymalnie cztery kopie, bez N=13 i bez otwierania task335.

Nie ponawiaj PHYS-001–007. Nie wykonuj PHYS-008/GEN-008 samodzielnie.
Nie zmieniaj A7, protokołu 333 ani statusu bariery. Nie kupuj hardware.

Dostarcz:
HANDOFF.md,
CLAIM_REGISTER.json z osobnym execution_status i scientific_verdict,
INPUTS.sha256 / MANIFEST.sha256,
minimalny replay,
NONCONCLUSIONS.md,
kontrprzykłady oraz tabelę INPUT/DERIVED/FITTED/TESTED.

Pokaż dokładnie, jakie nowe operatory generuje eliminacja kopii,
czy mogą być kontrolowanie ograniczone i czy efekt wyróżnia FIN
wobec tych samych kontrmodeli.

Po GEN-003 zatrzymaj się i oddaj checkpoint architektowi.
Jeżeli potrzebujesz nowej zasady fizycznej, nie wymyślaj jej sam.
~~~

## 10. Późniejsza droga i warunki STOP

Dopiero po wcześniejszych bramkach sensowna jest kolejność:

source/primitive law
→ autonomia i kompozycja
→ operacyjna lokalność
→ scaling i continuum
→ zachowane pola/propagacja
→ stabilne wzorce materii
→ osobne mosty quantum/gravity
→ niezależne testy fizyczne.

Jeżeli parent już zakłada przestrzeń, zegar, quantum dynamics albo gravity,
nie liczyć ich jako wygenerowanych wyników.

Zatrzymać daną interpretację, gdy:
- A7 jest jedynie programowanym wejściem bez niezależnego source;
- effective action jest dopasowywany od nowa na każdym poziomie;
- locality pojawia się tylko po ręcznym wpisaniu grafu;
- chronione wielkości wynikają z ręcznie ograniczonej gate class bez source;
- pamięć jest porównywalna z następną skalą, lecz pomijana;
- kontrmodel daje te same obserwacje w granicach błędu;
- poprawny wynik inżynierski jest używany jako dowód fundamentalności.

Najbliższy uczciwy cel nie brzmi „wygenerować rzeczywistość w symulacji”.
Brzmi: ustalić, czy z jednego jawnego FIN można wyprowadzać kolejne
opisy fizyczne, z kontrolą błędu i niezależną treścią predykcyjną.
