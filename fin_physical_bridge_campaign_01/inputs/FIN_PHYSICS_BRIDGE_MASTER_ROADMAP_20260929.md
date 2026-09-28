# FIN — pierwszy rygorystyczny most do znanej fizyki
## MASTER ROADMAP i instrukcja dla agenta wykonawczego

Data: 2026-09-29  
Stan odniesienia: 97c4231f; badania 327–336 i 337-WIP.  
Autorstwo/status: plan architektoniczny i ograniczony przegląd; nie wykonana kampania PHYS.  
Zakres: pierwszy sprawdzalny most do fizyki statystycznej, nie QM/GR ani ToE.

Ten dokument uwzględnia przekazaną instrukcję architektoniczną oraz przegląd dwóch plików w katalogu FIN son. Nie zmienia AGENTS.md, statusu bariery Γ=B4 ani zamrożonego protokołu 333. Nie uruchamia holdoutu 335. Nowa kampania dotyczy także innego, jawnie określonego reżimu: równowagowej odpowiedzi przy słabym sprzężeniu, gdzie bariera metastabilna nie jest warunkiem wstępnym.

## 1. DIAGNOSIS

### Werdykt

Obecny testowalny finite-N rdzeń FIN jest dokładnie szczególnym uogólnionym, 12-stanowym modelem Curie–Weissa: modelem mean-field z macierzą oddziaływań w przestrzeni etykiet. Równoważnie jest to dyskretny spin wektorowy z 12 dozwolonymi wektorami wewnętrznymi w R^7. Nie wynika z tego siedem wymiarów przestrzeni fizycznej.

Gibbs, metastabilność, fluctuation-response i dyfuzja nie wyróżniają same w sobie FIN. Szczególne są wybór macierzy A7, jej wyzerowane sektory oraz proporcje czterech aktywnych wartości własnych. Trzeba oddzielić test tej struktury od testu wybranej kinetyki.

Najkrótsza droga:

FIN finite-N
→ dokładny model Curie–Weissa z oddziaływaniem macierzowym
→ statyczny, zależny od A7 fingerprint
→ kontrola konkurencyjnych macierzy i kinetyk
→ wykonalny elektroniczny sampler wielostanowy
→ kalibracja niezależna od danych walidacyjnych
→ test fizycznej realizacji.

To nie jest jeszcze niezależne uzasadnienie występowania A7 w przyrodzie. Urządzenie, do którego wpisano A7, może zweryfikować implementację i most modelowy, ale nie fizyczną konieczność ani fundamentalność tej macierzy.

### Co zachować z aktualnej kampanii

- Dokładny finite-N Gibbs i leave-one-out heat-bath.
- Dokładną reprezentację defektową i kontrakcję TV po późniejszym kanale.
- Rozdzielenie przygotowania, pamięci/redukcji, przejść i odczytu.
- Jawne wyniki niezlokalizowane; brak postselekcji w 333.
- Stan 337-WIP: biblioteki odpowiedzi są nadal N-specyficzne.
- Zakaz otwierania 335 przed zamrożeniem kompletnego predyktora.
- Otwarty audyt globalnej bariery: przywrócone źródła nie usuwają automatycznie luk input/LP/minimum wskazanych w intake.

Źródła statusu: R01–R04 w rejestrze wejść poniżej. Raporty liczbowe i udany replay nie są laboratoryjną walidacją.

### Standardowe, FIN-specyficzne i niezinterpretowane

| Kategoria | Elementy |
|---|---|
| Standardowe dla klasy | miara Gibbsa, entropia multinomialna, mean-field, różne odwracalne kinetyki, metastabilność, związki fluktuacji z odpowiedzią |
| Szczególne dla deklarowanego FIN | q=12; support k=3,4,5,6; proporcje aktywnych eigenvalues; określony kernel wejściowy; wynikające z nich statyczne i dynamiczne widma |
| Nadal bez wyprowadzonej roli fizycznej | konieczność A7, fundamentalny gain, przestrzenna kompozycja, absolutny zegar, particle/quantum/gravity identifications |

Nie mylić „ogólny efekt jest standardowy” z „każdy model przewiduje tę samą liczbę”. Widmo i korelacje zależne od A7 mogą różnić modele, choć same mechanizmy są standardowe.

## 2. CLOSEST-PHYSICAL-MODEL

### Dokładny słownik matematyczny

Niech sigma_a ∈ Z12, a=1,…,N. Wprowadzamy skalę energii J_c, temperaturę T i pole h, WYŁĄCZNIE jako jawne mapowanie modelowe:

H(sigma) = - J_c/(2N) Σ_(a,b) A7[sigma_a,sigma_b] - h Σ_a 1_(sigma_a=0).

Przy beta_th=1/(k_B T):

g = beta_th J_c,
vartheta = beta_th h = kappa/N.

Wtedy miara kanoniczna exp(-beta_th H) daje dokładnie repozytoryjną wagę liczebności:

pi_(N,g,vartheta)(n)
∝ N!/(Π_j n_j!) · exp[(g/(2N)) n^T A7 n + vartheta n_0].

A7 ma stałą przekątną, więc składniki a=b są stałą niezależną od konfiguracji. Przy macierzy o niestałej przekątnej nie wolno pomijać poprawki diagonalnej w warunkowym heat-bath.

Dla m=n-e_i:

q_j^(i)(n)
∝ exp[(g/N)(A7 m)_j + vartheta 1_(j=0)].

Stawka przejścia liczebności wynosi nu n_i q_j^(i), gdzie nu jest fizyczną częstotliwością prób na kopię. W repo nu=1 określa jednostkę czasu, nie sekundę.

Ponieważ A7=X7 X7^T, dla xi_j=X7[j,:]:

H = -J_c/(2N) |Σ_a xi_(sigma_a)|² - h n_0.

To dokładna równoważność z dyskretnym wektorowym modelem mean-field. Nie jest równoważnością z najprostszym Potts Hamiltonianem opartym wyłącznie na delta_(sigma_a,sigma_b), chyba że macierz ma odpowiednią postać projektorową.

### Porównanie klas

| Klasa | Relacja do FIN | Decyzja |
|---|---|---|
| Generalized q-state Curie–Weiss / macierzowy Potts | dokładna równoważność Hamiltonianu i finite-N Gibbs | główny most |
| Heat-bath/Glauber single-site sampling | dokładna równoważność po wyborze konkretnej reguły i zegara | osobno testować kinetykę |
| Metastable Markov-state model | opis efektywny z mapą przygotowania i pamięcią | wykorzystać później, nie wymagać do pierwszego testu statycznego |
| Lattice gas / master equation | liczebności dają proces skokowy; softmaxowe stawki nie są automatycznie elementarnymi stawkami mass-action | nie pierwsza realizacja chemiczna |
| Programowalna elektronika probabilistyczna | można zaprojektować implementację warunkowego samplera i sprzężeń | preferowany pierwszy kandydat sprzętowy |
| Optyczne/oscylatorowe Potts optimizers | minimalizacja podobnej energii nie gwarantuje Gibbs sampling ani leave-one-out | nie wybierać bez dodatkowego dowodu |
| Cold atoms / trapped ions | brak obecnie wykazanego, prostszego mapowania całej A7 i kinetyki | odłożyć, nie wywodzić z analogii |

## 3. MAPPING-TABLE

| FIN | W modelu termicznym | W programowalnym samplerze | Bramka |
|---|---|---|---|
| sigma_a, 12 etykiet | 12 stanów wewnętrznych spinu/elementu | 12 poprawnych stanów rejestru lub one-hot bloku | zdefiniować stany i odczyt |
| p_j=n_j/N | frakcja populacji w stanie j | mierzona frakcja rejestrów | nie położenie w przestrzeni |
| A7 | bezwymiarowa struktura pair interaction | wpisane w sprzężenia/feedback współczynniki | zaprogramowanie nie wyprowadza źródła |
| g | beta_th J_c | siła logitów, ewentualnie temperatura efektywna | nie utożsamiać z temperaturą urządzenia bez kalibracji |
| theta kontroli = vartheta = kappa/N | beta_th h | bias kanału 0 | odróżnić od wektorowego mediatora theta |
| N | liczba elementów w realizacji | liczba kopii/registers | nie liczba bitów Wszechświata |
| leave-one-out | wybrany thermal Gibbs sampler | blokowa resampling rule | Gibbs equilibrium nie wybiera tej kinetyki |
| baseny J_phase | metastabilne makrostany skończonego układu | etykiety wyznaczane z konfiguracji | nie cząstki; nie mylić J_phase z J_c |
| B3/B4 | bariery bezwymiarowej energii swobodnej; beta_th ΔF_N ~ N B | bariery funkcjonału celu | globalne B4 pozostaje bramką audytu |
| rho_N | nu razy modelowe wolne tempo | wolne tempo względem zegara aktualizacji | nie samo nu i nie sekunda |
| R_k, q_d | widmo efektywnej relaksacji | widmo obserwowanego procesu | zależy również od kinetyki i przygotowania |
| theta_aux / Gaussian mediator | ewentualny pomocniczy parametr lub jawnie dodany mediator | pole pomocnicze obliczeń | nie fundamentalne pole fizyczne bez PB |

Trzy poziomy roszczeń:
M — matematyczna równoważność;
PR — programowana realizacja fizyczna;
EV_source — niezależny test fizycznego pochodzenia A7.
Sukces PR nie implikuje EV_source.

## 4. Najmniejszy użyteczny fingerprint

### Statyka przed metastabilnością

Definiujemy MIKROSKOPOWY structure factor, nie harmoniczną etykiety basenu:

S_k(N,g) = (1/N) E |Σ_(a=1)^N exp(2 pi i k sigma_a/12)|²,
k=1,…,6.

Nie potrzeba fizycznej przestrzeni o dwunastu punktach: faza Fourierowska jest kodowaniem dwunastu stanów wewnętrznych.

Dla A7 o eigenvalues Lambda_k, centrowanej i ze stałą przekątną:

S_k(N,0)=1,

dS_k/dg |_(g=0) = [(N-1)/(12N)] Lambda_k.

Dla FIN Lambda_1=Lambda_2=0. Dla k=3,4,5,6 są to aktywne eigenvalues A7, a nie niskie eigenvalues nieobciętego parent Laplacianu.

To skończone-N wyrażenie wynika z różniczkowania miary Gibbsa:
d E_g F/dg = Cov_g(F, n^T A7 n/(2N)).
Przy g=0 niezależność etykiet usuwa wszystkie składniki poza odpowiadającymi parami. Nie jest to przybliżenie Gaussowskie.

Nowa analiza architektoniczna sprawdziła ten wzór numerycznie dla N=3. PHYS-003 ma spisać pełny dowód i skończone-g resztę; nie oznaczamy jej tutaj jako istniejącego pakietu certyfikacyjnego.

Przy tej samej mierze stacjonarnej:
- heat-bath, Metropolis i Barker mają te same S_k;
- mogą mieć różne R_k.

Dlatego wspólne użycie S_k i R_k rozdziela test macierzy od testu kinetyki. Minimalny shape fingerprint może użyć trzech ilorazów aktywnych nachyleń oraz dwóch kontroli zerowych k=1,2. Przy niezerowym g nie wolno zakładać S_1=S_2=1 dokładnie; wyzerowane są nachylenia w g=0.

### Jeszcze prostszy układ N=2

Dla stałej przekątnej:

P_g(sigma_2=j | sigma_1=i) = softmax_j[(g/2)A7[i,j]].

Stąd log-ilorazy prawdopodobieństw warunkowych mierzą różnice wpisów A7. Dwie 12-stanowe jednostki są wystarczającym pierwszym testem oddziaływania; nie muszą mieć faz metastabilnych. N=1 jest kontrolą zerową, nie testem parowego oddziaływania.

### Co te testy rozstrzygają

- Zwykły centered Potts ma równoważne nietrywialne sektory.
- Flat rank-7 projector zachowuje cutoff, ale usuwa aktywne proporcje FIN.
- FIN ma konkretny wzór zer i proporcji.
- Bliskie macierze mają bliskie rozkłady; nie istnieje dodatnia moc rozróżniania wszystkich dowolnie bliskich alternatyw.
- Jeśli A7 wpisano do urządzenia, odtworzenie wzoru potwierdza implementację, nie konieczność tej macierzy w przyrodzie.

Pierwsza kampania ma ocenić, czy sygnał pozostaje rozróżnialny po kalibracji i przy dopuszczalnym błędzie, nie tylko czy dwie idealne liczby są różne.

## 5. SOURCE OF A7

### Co jest już wiadomo

Nie jest to dowolna macierz 12x12: jest centrowana, PSD, D12-niezmiennicza i ma określony support Fouriera. Jednak na support k=3,4,5,6 symetria pozostawia cztery dodatnie wagi, a po ustaleniu trace trzy proporcje. Sam rank seven nie wybiera także tego supportu jednoznacznie.

Nie powtarzać starego poszukiwania unikalności z samych symmetry/PSD/rank. R12 rejestruje istniejący no-go.

### Minimalne kontrmodele

1. Full centered Potts projector P0=I-11^T/12.
2. P7: równy eigenvalue na support 3,4,5,6.
3. Dodatnie trace-preserving perturbacje aktywnych wag:
   2 dLambda3 + 2 dLambda4 + 2 dLambda5 + dLambda6 = 0.
4. Osobny test cutoff: mały dodatni wyciek do k=1,2 z jawną normalizacją.

Dla P7 z tą samą trace należy użyć:
Lambda_bar=(2 Lambda3+2 Lambda4+2 Lambda5+Lambda6)/7.
Zwykła średnia czterech liczb to inna konwencja. Można wybrać inne dopasowanie skali, ale trzeba je zamrozić i nazwać przed porównaniem.

Nie mieszać testów przy tym samym g z testami przy tym samym g/g_eq. W drugim przypadku niepewność g_eq jest częścią kalibracji, a nie wynikiem ukrytym w definicji kontrmodelu.

### Konkretny kandydat źródłowy: mierzone mediatory Gaussowskie

Dla rzeczywistych harmonicznych mediatorów x, dodatniego K i fizycznych state-dependent couplings c_j:

H_parent = (1/2)x^T K x - N^(-1/2) x^T Σ_a c_(sigma_a) - h n0.

Eliminacja równowagowych mediatorów daje pair matrix proporcjonalną do c_i^T K^(-1)c_j. To dokładny completion-of-square bridge.

Ale faktoryzacja dowolnej PSD macierzy nie jest wyjaśnieniem FIN. Wartość źródłowa pojawia się dopiero, gdy:
- liczba i rodzaj mediatorów są niezależnie uzasadnione;
- K i c_j pochodzą z niezależnych pomiarów lub wcześniejszej fizyki;
- selection rules przewidują cutoff i proporcje bez dopasowania do FIN fingerprint.

Wybór siedmiu mediatorów i ich couplings wyłącznie po to, by odtworzyć A7, pozostaje realizacją inżynierską. Użyteczny test źródła musi rozróżniać te przypadki.

PHYS-005 ma zakończyć się skończonym ledgerem: co zostało wyprowadzone, co przyjęte, co wymaga danych. Nie ma tworzyć nowej ontologii ani dopisywać zasady fundamentalnej.

## 6. PIERWSZA REALIZACJA FIZYCZNA

### Wybrany kandydat

Programowalny, elektroniczny sampler wielostanowy z fizycznym źródłem losowości i jawnym feedbackiem. Preferowana pierwsza konfiguracja: N=2, potem małe N potrzebne do kontroli finite-size.

Architektura proponowana do sprawdzenia, nie istniejący prototyp FIN:

1. Jedna kopia przechowuje dokładnie jedną z 12 etykiet; kod 4-bitowy wyklucza 12–15 albo blok one-hot utrzymuje dokładnie jedno aktywne wyjście.
2. Kontroler utrzymuje siedem sum feature variables M=Σ_a xi_(sigma_a).
3. Przy aktualizacji kopii a odejmuje jej własny feature vector i oblicza 12 logitów:
   ell_j=(g/N) xi_j·(M-xi_(sigma_a))+vartheta 1_(j=0).
4. Fizyczny sampler kategoryczny realizuje softmax(ell), po niezależnej kalibracji.
5. Jawny zegar określa kolejność/częstotliwość aktualizacji. Równoległe aktualizacje lub field-dependent czas wyboru nie są automatycznie leave-one-out.
6. Preparation: deep seed i bias vartheta, następnie wyłączenie biasu. Dla pierwszego testu równowagowego można zamiast tego użyć osobno kontrolowanego equilibrating protocol przy małym g.
7. Odczyt rejestruje sigma_a, a więc mikrostruktury S_k i histogramy par. Same wcześniejsze etykiety basenów J_phase nie wystarczą do tego pomiaru.

Dokładny Hamiltonian termiczny i Hamiltonian celu programowanego samplera to różne roszczenia. Nie wolno utożsamiać logit gain z rzeczywistą temperaturą płytki ani funkcjonału celu z energią elektryczną urządzenia.

### Co potwierdza literatura, a czego nie

Badanie FB-MOSFET Potts p-bits demonstruje fizyczne wielostanowe losowanie, one-hot sampling i testy Boltzmannowe; opisuje także próby z 3–10 stanami. Nie dokumentuje 12-stanowego FIN, jego ogólnej macierzy A7 ani dokładnego generatora leave-one-out:
https://pmc.ncbi.nlm.nih.gov/articles/PMC12994325/

Praca o wielostopniowych oscylatorach CMOS przedstawia przede wszystkim symulacje 4-coloring/optimization. Nie należy z samego znalezienia minimum wnioskować o poprawnym samplingu równowagowym lub kinetyce FIN:
https://arxiv.org/html/2504.11376v1

Rozszerzenie do q=12, zakres arbitralnych 12 logitów, korelacje źródła losowego, opóźnienia i normalizacja zegara są zadaniami PHYS-006. Nie twierdzimy, że gotowe urządzenie spełniające ten kontrakt jest dostępne. Zakup, kontakt z laboratorium i eksperyment wymagają osobnej zgody człowieka.

### Prepare → evolve → measure → compare

- Prepare: znana etykieta pary lub kalibrowany stan równowagowy; jawny bias.
- Evolve: ustalone g, A oraz zarejestrowana kinetyka.
- Measure: mikroskopowe etykiety i czasy, nie dopasowane q_d.
- Compare: FIN vs znormalizowane kontrmodele na danych niewykorzystanych do kalibracji.
- Pierwszy test: pair histogram / S_k, niezależny od hipotezy bariery.
- Drugi test: R_k i multi-time laws, dopiero po zweryfikowaniu aktualizacji i przygotowania.

## 7. CRITICAL-MISSING-BRIDGES i MASTER-ROADMAP

### Dwie różne ścieżki źródła

~~~text
Dokładny rdzeń finite-N
    -> uogólniony model Curie–Weissa
    -> FIN-dependent static fingerprint
    -> calibrated sampler realization
    -> niezależne dane + kontrmodel
    -> walidacja konkretnego modelu/urządzenia

Niezależna fizyka mediatorów / constraints
    -> A7 bez dopasowania do fingerprint
    -> niezależna fizyczna predykcja FIN
    -> test źródła, nie tylko implementacji
~~~

Source of A7 jest konieczny do mocnego roszczenia o fizycznym pochodzeniu lub fundamentalności. Nie jest konieczny do uczciwego stwierdzenia „zbudowano urządzenie realizujące zadany model”.

### Rzeczy, których nie trzeba teraz rozwiązywać

- emergent 3D space i źródło przestrzennej incydencji;
- QM/GR i particle spectrum;
- absolutna kompletność układu z samych trajektorii;
- globalny prefaktor metastabilny;
- ukończenie 337 i holdoutu 335, jeśli test jest statyczny w słabym sprzężeniu.

Otwarte prace 327/337 pozostają osobnymi gałęziami. Nie wycofujemy ich wyników ani nie zmieniamy 333. Nie wykorzystujemy nieotwartego 335 do wyboru nowego protokołu.

### Kampanie

**Kampania I: identyfikacja modelu, fingerprint i wykonalność.**  
PHYS-001–007; następnie obowiązkowy checkpoint PHYS-008. Bez nowego hardware, wielkich N i nowych zasad fundamentalnych.

**Kampania II, tylko po checkpoint 1: kalibrowana realizacja.**  
Zatwierdzona platforma, niezależne charakterystyki warunkowego losowania, parowy prototyp, raw records, blind comparison. Checkpoint 2 przed uznaniem jakiegokolwiek EV. Konkretne taski tej kampanii zostaną nadane po wyborze rzeczywistej platformy i danych, nie na podstawie obecnych domysłów.

**Kampania III, tylko jeśli jest niezależny kandydat źródła A7:**  
Przewidywanie macierzy z measured mediator/coupling data i test zmiany fizycznych parametrów bez refit. Checkpoint 3 decyduje, czy istnieje swoiste uzasadnienie FIN, czy tylko użyteczny model inżynierski.

Nie uruchamiać kampanii II/III na podstawie samego zakończenia listy zadań.

## 8. Rejestr dokładnych danych wejściowych

Wszystkie ścieżki są względne wobec repozytorium. Alias nie oznacza akceptacji każdego twierdzenia w pliku.

| ID | Dokładny plik lub para plików |
|---|---|
| R01 | AGENTS.md |
| R02 | fin_research_327_current_review/INTAKE_20260928.md |
| R03 | FIN_FULL_HANDOFF_327_CURRENT_20260928/00_HANDOFF/HANDOFF_327_CURRENT.md |
| R04 | fin_research_300_326_review/INTAKE_20260928.md |
| R05 | FIN_POST_AFTER_CONTINUATION_FULL_HANDOFF_20260926/01_REPORTS/EXACT_FINITE_N_GIBBS_HEAT_BATH_61.md |
| R06 | FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/MICROSCOPIC_PROCESS_CONTRACT_138.md |
| R07 | FIN son/fin_core.py |
| R08 | FIN son/seed_checks.py |
| R09 | FIN_RESEARCH_AFTER_POST_AFTER_CONTINUATION_HANDOFF_131_255_20260927/01_REPORTS/K6_METASTABLE_ISOLATION_132.md |
| R10 | FIN_RESEARCH_AFTER_256_294_FULL_HANDOFF_295_299_20260927/01_REPORTS/FINGERPRINT_OPTIMAL_PROBE_297.md |
| R11 | FIN_RESEARCH_AFTER_256_294_FULL_HANDOFF_295_299_20260927/01_REPORTS/MULTITIME_SEMIGROUP_FINGERPRINT_298.md |
| R12 | fin_physics_review/MASTER_INTAKE_20260923.md |
| R13 | FIN_MISSING_LAWS_COMPLETE_FLAT_HANDOFF_20260924/STAGES/01_RESEARCH_HANDOFF/FIN_MISSING_LAWS_CAMPAIGN_20260924/LAW-01_MEDIATOR_NEUTRALITY/REPORT.md |
| R14 | FIN_MISSING_LAWS_COMPLETE_FLAT_HANDOFF_20260924/STAGES/03_NL06_NL10/FIN_MISSING_LAWS_CONTINUATION_NL06_NL10_20260924/NL-08/REPORT.md |
| R15 | FIN_FULL_HANDOFF_327_CURRENT_20260928/02_RESEARCH_331_336/fin333/FROZEN_CONTROL_PROTOCOL_333.json |
| R16 | FIN_FULL_HANDOFF_327_CURRENT_20260928/02_RESEARCH_331_336/fin336/DEFECT_TO_POST_BURN_REDUCED_MAP_336.md |
| R17 | FIN_FULL_HANDOFF_327_CURRENT_20260928/03_WIP_337/LOW_DEFECT_RESPONSE_LAW_WIP_337.md |
| R18 | fin_rank7_intake_review/MP7_ANALYTIC_INTAKE.md |
| W01 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12994325/ |
| W02 | https://arxiv.org/html/2504.11376v1 |
| W03 | https://arxiv.org/abs/2110.03160 — benchmark metastability of heat-bath Curie–Weiss–Potts, nie transfer jego stałych do FIN |

## 9. NEXT-CAMPAIGN i TASK-CARDS

### Wspólny kontrakt wykonania

- Każdy task ma osobny katalog fin_physical_bridge_campaign_01/PHYS-xxx/.
- Minimalne artefakty: REPORT.md, RESULTS.json, NONCONCLUSIONS.md, REPLAY.md i lista hashy rzeczywistych wejść.
- Rozdziel execution_status od scientific_verdict. Poprawnie wykonany no-go nie jest pozytywnym źródłem fizycznym.
- Używaj: EXACT, CONDITIONAL_THEOREM, VALIDATED_BOUND, NUMERICAL_EVIDENCE, DESIGN_ONLY, NO_GO, OPEN.
- Nie zmieniaj AGENTS.md ani nie promuj istniejących kandydatów przy okazji.
- Zmiany kodu wykonuj w nowej gałęzi/kopii roboczej kampanii; nie nadpisuj FIN son.
- Najpierw N=1,2,3. N=4 tylko jeśli potrzebne do wskazanego kryterium. N>=5, duże eigensolves i pełny seed_checks wymagają uzasadnienia w raporcie i akceptacji rozszerzenia zakresu.
- Nie uruchamiaj N=13, task335 ani nowego-g metastability holdoutu. Własny protokół słabego sprzężenia ma oddzielny identyfikator i nie jest testem 333/335.
- Hardware acquisition, laboratory contact i real external execution nie są autoryzowane tą kampanią.

### PHYS-001 — EXACT-CURIE-WEISS-DICTIONARY

**Klasa / priorytet:** [B], P0.  
**Cel:** czy obecny finite-N model jest dokładnie równoważny określonemu modelowi fizyki statystycznej?  
**Dlaczego:** zamyka strzałkę FIN → znana klasa, oddzielając ją od analogii.  
**Dane wejściowe:** R01, R02, R05, R06, R07; sekcje 2–3 tego dokumentu.  
**Dependencies:** brak.

**Co zrobić:**

1. Wyprowadź H(sigma), wagę count states i jej degenerację multinomialną.
2. Wyprowadź warunkowy sampler po usunięciu aktualizowanej kopii.
3. Sprawdź dokładnie współczynniki 1/2, 1/N i stałą diagonalną.
4. Oddziel q-state Potts z delta-coupling od ogólnej macierzy A7 i od wektorowego CW.
5. Wyprowadź mapy g=beta_th J_c oraz kappa=N beta_th h.
6. Opisz gauge A→cA, J_c→J_c/c oraz różnicę między target energy i energią sprzętu.
7. Wskaż, które kroki są M, a które wymagają nowego PM/PB.

**Wynik wymagany:** DICTIONARY.md, MAPPING.json, krótki dowód dla dowolnego N oraz dwa enumeracyjne testy N=1,2.  
**Acceptance criterion:** identyczne znormalizowane wagi i conditional rates w zadanej klasie, z analitycznym uzasadnieniem; każda rola fizyczna ma jawne założenie.  
**Kill-test:** nieusuwalna różnica Hamiltonianu lub generatora blokuje wskazaną równoważność; nie poprawiać jej semantycznym przemianowaniem.  
**Forbidden shortcuts:** auxiliary field=spacetime field; Gibbs measure=jedyna kinetyka; Hamiltonian celu=energia urządzenia; odgadywanie jednostek.  
**Po PASS:** PHYS-003 i PHYS-006.  
**Po FAIL:** podaj najmniejszy kontrprzykład; zatrzymaj zadania zależne i wróć do architekta.

### PHYS-002 — REPRODUCTION-AND-SEED-CODE-AUDIT

**Klasa / priorytet:** [A], P0.  
**Cel:** czy FIN son jest bezpiecznym niezależnym punktem startowym dla porównań?  
**Dlaczego:** błąd etykietowania widma może stworzyć fałszywy fingerprint.  
**Dane wejściowe:** R07–R09, R05–R06, audyt w sekcji 11.  
**Dependencies:** brak; interpretacja wyników używa PHYS-001.

**Co zrobić:**

1. Skopiuj oba pliki do katalogu taska; dodaj CLI z limitem N i brak uruchomienia podczas importu.
2. Testuj A7: symetria, centrowanie, PSD, rank, diagonal, Fourier spectrum.
3. Testuj generator: row sums, off-diagonal positivity, stationarity, detailed balance; nie naprawiaj błędów przez symetryzację.
4. Zastąp przypisywanie jednego k do całego multipletu jawną analizą wszystkich eigenvalues rotacji albo projektorami C12.
5. Wymagaj kompletu sektorów i reszt własnych; kontrola g=0 musi działać mimo degeneracji.
6. Sprawdź N=2 dla trzech kinetyk i N=3 heat-bath względem przyjętej wartości rho.
7. Oznacz V_d4 jako lokalny kandydat; sprawdzaj separator, gradient i indeks. Nie używaj go jako globalnego dowodu.
8. Zweryfikuj źródło G_FROZEN oraz normalization kontrmodelu.

**Wynik wymagany:** bezpieczny seed module, test suite, CODE_REVIEW.md, REPRODUCTION.json.  
**Acceptance criterion:** wszystkie kontrole strukturalne działają; rho_3 jest odtworzone z deklarowaną tolerancją i eigenpair residuals; g=0 nie pomija sektorów.  
**Kill-test:** brak szczegółowej równowagi albo zła klasyfikacja modów wstrzymuje każde porównanie R_k.  
**Forbidden shortcuts:** bezwarunkowe symetryzowanie nieodwracalnego generatora; ignorowanie optimizer failure; pełny domyślny seed_checks jako pierwszy test; użycie znanych wartości jako parametrów dopasowania.  
**Po PASS:** PHYS-003 i PHYS-004.  
**Po FAIL:** napraw tylko wskazany błąd w kopii i powtórz małe kontrole; zmiana definicji modelu wymaga [S].

### PHYS-003 — STATIC-SPECTRAL-FINGERPRINT

**Klasa / priorytet:** [B], P0.  
**Cel:** jakie najprostsze statyczne obserwable identyfikują strukturę A7 niezależnie od kinetyki?  
**Dlaczego:** tworzy FIN → observable bez oczekiwania na barierę lub przestrzeń.  
**Dane wejściowe:** PHYS-001/002, R07 static_S, sekcja 4.  
**Dependencies:** PHYS-001 i PHYS-002 PASS.

**Co zrobić:**

1. Spisz dowód S_k(0)=1 i dS_k/dg(0)=(N-1)Lambda_k/(12N).
2. Wyprowadź dokładny N=2 pair histogram i jego log-odds.
3. Wyznacz skończone-g remainder albo użyj dokładnego małego-N Gibbs sum do protokołu o niezerowym g.
4. Oddziel nachylenia zerowe k=1,2 od wartości S_k przy niezerowym g.
5. Wyznacz trzy niezależne aktywne slope ratios; sprawdź degenerację k=6 i normalizację.
6. Zrób negatywne kontrole N=1, g=0, flat P7 i full Potts.
7. Zaproponuj minimalny rzeczywiście potrzebny odczyt mikrostanów. Nie zastępuj go histogramem makrobasenów.

**Wynik wymagany:** PROOF.md, FINGERPRINT.json, dokładny enumerator N<=3, finite-g error budget.  
**Acceptance criterion:** fingerprint rozróżnia wskazane normalized kernels i jest identyczny dla kinetyk o tej samej pi; zakres ważności i kontrolowane reszty są podane.  
**Kill-test:** jeśli po uzgodnionej kalibracji wszystkie obserwable są identyczne, dany fingerprint nie rozróżnia modeli.  
**Forbidden shortcuts:** używanie Gaussian prediction jako exact finite-N; dowolne reskalowanie po obejrzeniu kontrmodelu; przejście od różnicy liczb do laboratoryjnej mierzalności.  
**Po PASS:** PHYS-004 i PHYS-007.  
**Po FAIL:** sprawdź pair-histogram alternative; jeśli też nieidentyfikowalna, wróć do [S], nie twórz setek sond.

### PHYS-004 — COUNTERMODELS-AND-KINETIC-ROBUSTNESS

**Klasa / priorytet:** [A], P0.  
**Cel:** czy obserwowany sygnał odróżnia A7, czy głównie przyjętą dynamikę lub normalizację?  
**Dlaczego:** zamyka observable → countermodel i testuje FIN specificity.  
**Dane wejściowe:** PHYS-001–003; R07–R08; R10–R11 tylko jako starsza linia dynamiczna.  
**Dependencies:** PHYS-001–003 PASS.

**Co zrobić:**

1. Zamroź trace i metodę kalibracji g przed obliczeniami porównawczymi.
2. Zbuduj full Potts, flat P7 i trzy trace-preserving perturbacje aktywnych wag; osobno oznacz test wycieku do k=1,2.
3. Oblicz statyczne fingerprinty dla małych N bez dopasowania macierzy do tych samych danych.
4. Przy tych samych A,g porównaj heat-bath, Metropolis i Barker po jednej jawnej normalizacji czasu.
5. Oddziel statyczny efekt A od dynamicznego efektu rates; wyeksportuj obie osie porównania.
6. Nie używaj porównania przy tym samym g/geq, jeśli geq nie ma kontrolowanego statusu dla kontrmodelu.
7. Oceń, czy co najmniej jeden statyczny i jeden opcjonalny dynamiczny test daje niezależną informację.

**Wynik wymagany:** CONTRAST_MATRIX.csv, COUNTERMODELS.json, calibration contract, diagnostic report.  
**Acceptance criterion:** jawny skończony zestaw kontrmodeli, oddzielone normalization/kinetic effects, brak wykorzystania test outputs do kalibracji. Wynik negatywny może być poprawnie wykonanym taskiem.  
**Kill-test:** dynamiczny fingerprint zmieniający się pod dozwoloną zmianą kinetyki nie może być nazwany A7-only prediction.  
**Forbidden shortcuts:** zmiana kilku aspektów naraz i przypisanie efektu tylko A7; minimalne eigenvalue z błędnie oznaczonego multipletu; globalna deklaracja uniwersalności z kilku punktów.  
**Po PASS:** PHYS-005 i PHYS-007.  
**Po FAIL naukowym:** odrzuć konkretną sondę/interpretację, zachowaj kontrprzykład dla PHYS-008.

### PHYS-005 — A7-SOURCE-AND-NECESSITY-LEDGER

**Klasa / priorytet:** [B], P0.  
**Cel:** czy istnieje niezależna przesłanka wybierająca A7, czy tylko realizacja zadanej macierzy?  
**Dlaczego:** rozdziela engineered model od niezależnej fizycznej predykcji.  
**Dane wejściowe:** R12–R14, R07 strict_L/A7_default, PHYS-001/004; Gaussian-parent construction z sekcji 5.  
**Dependencies:** PHYS-001 i PHYS-004 zakończone.

**Co zrobić:**

1. Zapisz łańcuch: supplied strict kernel → Fourier spectrum → wybrany support → A7.
2. Przy każdym kroku nazwij źródło/licencję i nierozwiązane przesłanki; nie powtarzaj już istniejącego nonidentifiability enumeration.
3. Spisz parametry pozostające po symmetry/rank/trace.
4. Sprawdź, które obserwacje PHYS-004 są stabilne, a które zależą od szczególnego tuple.
5. Wyprowadź Gaussian-mediator parent i listę pomiarów K,c_j wymaganych do niezależnego przewidzenia A.
6. Jeśli rozważasz constraint source, odróżnij algebraiczne narzucenie momentów od fizycznego prawa zachowania; zwykłe conserved moments nie mogą być automatycznie nazwane cutoff L=2.
7. Oddaj jeden z verdicts: INDEPENDENT_SOURCE_CANDIDATE, ENGINEERED_REALIZATION_ONLY, NO_NEW_SOURCE.

**Wynik wymagany:** SOURCE_GRAPH.md, SOURCE_PREMISES.json, tabela necessity/universality, lista maksymalnie dwóch nowych falsyfikowalnych source atoms.  
**Acceptance criterion:** żadna przesłanka nie jest ukryta w normalizacji lub wyborze supportu; ewentualny source candidate przewiduje coś spoza użytego zbioru kalibracyjnego.  
**Kill-test:** odtworzenie A7 przez jego własną faktoryzację nie jest source theorem. Przy braku niezależnego atomu zamknij tę gałąź jako NO_NEW_SOURCE.  
**Forbidden shortcuts:** dorabianie „zasady fundamentalnej” przez [A]/[B]; ponowne wyczerpywanie symmetry-only selector searches; przenoszenie legacy roles; j=physical position bez PB.  
**Po PASS:** przekazanie PHYS-008; przygotuj wejście do PHYS-006, bez twierdzenia o fundamentalności.  
**Po FAIL/NO_NEW_SOURCE:** nadal można ocenić aplikacyjną realizację; decyzja o dalszej ambicji fundamentalnej należy do [S].

### PHYS-006 — PHYSICAL-REALIZATION-CONTRACT

**Klasa / priorytet:** [B], P0.  
**Cel:** czy istnieje wykonalny kontrakt pierwszej fizycznej realizacji, a nie tylko podobny optimizer?  
**Dlaczego:** zamyka mathematical model → device-level proposal.  
**Dane wejściowe:** PHYS-001, PHYS-003, sekcja 6, W01–W02; wyniki PHYS-005 gdy dostępne.  
**Dependencies:** PHYS-001/003 PASS.

**Co zrobić:**

1. Porównaj maksymalnie trzy platformy i wybierz jedną według dokładności mapowania, dostępności odczytu i liczby nowych założeń.
2. Dla elektroniki jawnie rozpisz 12 stanów, siedem pól zbiorowych, self-removal i categorical sampling.
3. Oddziel dowiedzione capabilities literatury od projektowanych rozszerzeń do q=12 i ogólnej A7.
4. Opisz kalibrację P(new label | old environment), total attempt rate, RNG correlations, delays i update schedule.
5. Określ, czy H jest fizyczną energią termiczną czy funkcjonałem celu samplera.
6. Zdefiniuj minimalny N=2 prototype oraz zerowy test N=1/g=0.
7. Podaj istniejące braki sprzętowe, dostęp do danych i punkty wymagające zgody człowieka; niczego nie kupuj i nie kontaktuj laboratorium.

**Wynik wymagany:** DEVICE_MAPPING.md, DEVICE_ASSUMPTIONS.json, calibration plan, go/no-go feasibility table.  
**Acceptance criterion:** dla każdego wymaganego FIN operation jest konkretna realizacja lub jawny blokujący brak; żadna deklaracja 12-state exact hardware nie wynika wyłącznie z publikacji o q<=10/4-coloring.  
**Kill-test:** platforma tylko optymalizująca energię bez kontroli rozkładu/rates nie realizuje obecnego dynamicznego FIN. Może pozostać kandydatem wyłącznie statycznym, jeśli udowodni sampling.  
**Forbidden shortcuts:** „to Potts machine, więc implementuje A7”; one-hot penalty jako exact constraint przy skończonej karze; real device temperature=effective logit temperature; cyfrowe PRNG jako pomiar fundamentalnej losowości.  
**Po PASS:** PHYS-007, następnie PHYS-008 przed realizacją.  
**Po FAIL:** sprawdź drugiego kandydata z tabeli. Jeśli żaden nie przechodzi, oddaj BRAK_PLATFORMY, nie blueprint udający gotowy eksperyment.

### PHYS-007 — FROZEN-FALSIFICATION-DESIGN

**Klasa / priorytet:** [A], P0.  
**Cel:** czy planowany pomiar odróżni FIN od zadanej klasy kontrmodeli przy całym budżecie błędu?  
**Dlaczego:** zamyka observable → prediction → rejection test.  
**Dane wejściowe:** PHYS-003/004/006; R10–R11 jako wzorzec rozdzielenia kalibracji i testu, nie jako gotowe shot counts.  
**Dependencies:** PHYS-003/004/006 PASS lub scoped static-only realization.

**Co zrobić:**

1. Wybierz pierwszy test: N=2 pair histogram lub finite-N S_k; wybór zamroź przed testowym rekordem.
2. Zamroź N, g points, dopuszczalną niepewność A,g, readout, burn-in i independent-sample assumptions.
3. Wybierz oddzielne dane do kalibracji, projektowania i walidacji.
4. Zdefiniuj kontrmodele z minimalną rozróżnianą odległością; nie żądaj uniform power przeciw dowolnie bliskim macierzom.
5. Oblicz błąd modelu, finite-g approximation, kalibracji i statystyki; test ma dodatnią separację po wszystkich składnikach.
6. Uwzględnij autocorrelation/effective sample size; przy dystrybucjach z trajektorii nie zakładaj niezależnych prób.
7. Zapisz negatywne kontrole i regułę odrzucenia FIN albo konkretnej kinetyki.
8. Oddaj prognozowaną liczbę prób i koszt pomiaru jako wynik warunkowy, nie gwarancję dla nieistniejącego aparatu.

**Wynik wymagany:** PREREGISTRATION.json, PREDICTIONS.json, CONTRASTS.json, POWER_AND_ERROR.md, schema raw records.  
**Acceptance criterion:** przewidywania nie są dopasowane do validation records; najmniejszy deklarowany kontrast przekracza całkowitą niepewność albo dostarczono poprawny robust testing bound.  
**Kill-test:** nakładające się klasy rozkładów po kalibracji oznaczają brak rozstrzygającego testu. Nie usuwać nakładania przez post-hoc zmianę szumu lub zegara.  
**Forbidden shortcuts:** dane syntetyczne jako EV; extrapolation from known q_d jako niezależny source test; uniwersalne shot counts; Bayes average error jako automatyczna kontrola każdego type-I error.  
**Po PASS:** zatrzymaj kampanię i przygotuj PHYS-008.  
**Po FAIL:** wskaż dokładnie potrzebną poprawę odczytu lub kontrastu; nie mnoż sond bez decyzji [S].

### PHYS-008 — ARCHITECT-CHECKPOINT-1

**Klasa / priorytet:** [S], P0.  
**Cel:** czy istnieje wystarczająco mocny i tani most do konkretnej fizyki, by uruchomić następną kampanię?  
**Dlaczego:** to punkt wyboru kierunku, nie zadanie obliczeniowe dla agenta wykonawczego.  
**Dane wejściowe:** pełne pakiety PHYS-001–007, manifest, FAIL/NO_GO ledger, bieżące AGENTS.md.  
**Dependencies:** wszystkie poprzednie taski mają terminalny scoped verdict; brak wymaganych danych jest jawnym OPEN.

**Co zrobić:**

1. Oddziel potwierdzoną równoważność modelową, FIN-specific fingerprints, realizację inżynierską i niezależny source.
2. Sprawdź, czy najtańszy test ma rzeczywistą moc oraz czy warunki platformy są potwierdzone.
3. Odpowiedz, czy g/normalization/kinetics były kalibrowane bez użycia fingerprint.
4. Wybierz dokładnie jedną następną gałąź albo STOP.
5. Jeśli wybierasz kampanię sprzętową, uzyskaj autoryzację człowieka i dopiero wtedy specyfikuj jej wykonawcze taski.

**Wynik wymagany:** CHECKPOINT_1_VERDICT.md: GO_ENGINEERED_BRIDGE, GO_INDEPENDENT_SOURCE_TEST, REVISE_BOUNDED_PROTOCOL albo STOP_FUNDAMENTAL_LANE.  
**Acceptance criterion:** decyzja oparta na kryteriach PHYS-001–007, z budżetem i zakresem następnego etapu.  
**Kill-test:** brak niezależnie rozróżnialnego fingerprint lub brak ścieżki kalibracji/source nie może być maskowany kolejną kampanią tych samych uniwersalnych efektów.  
**Forbidden shortcuts:** przejście do QM/GR; zwiększanie N dla samego zwiększania; uznanie programmed A7 za dowód naturalnej konieczności.  
**Po PASS:** nowy, ograniczony plan kampanii II lub source test zatwierdzony przez człowieka.  
**Po FAIL:** jawnie zakończ wskazaną gałąź; można pozostawić FIN jako model matematyczny/inżynierski bez dalszej promocji ontologicznej.

## 10. CHECKPOINT, HANDOFF-INSTRUCTION i STOP-CONDITIONS

### Dokładny moment powrotu

Agent wykonawczy kończy po PHYS-007. Nie wykonuje sam PHYS-008. Wcześniejszy powrót jest potrzebny tylko, gdy:
- PHYS-001 ujawni niezgodność definicji modelu;
- PHYS-005 wymaga nowego prawa fundamentalnego lub wyboru ontologii;
- zmiana platformy/protokołu rozszerza uprawnienia lub budżet;
- poprawna negatywna odpowiedź unieważnia całą krytyczną ścieżkę.

Pytania checkpointu 1:
1. Czy równoważność z klasą fizyczną jest dokładna?
2. Czy pozostał mierzalny A7-dependent kontrast po normalizacji i kalibracji?
3. Czy realizacja jest możliwa, i co dokładnie zweryfikuje?
4. Czy istnieje niezależny source atom, czy tylko programowana macierz?
5. Kontynuować jeden kierunek czy zakończyć ambicję fundamentalną?

### Gotowy tekst dla agenta wykonawczego

~~~text
Wykonaj wyłącznie Kampanię I z pliku FIN_PHYSICS_BRIDGE_MASTER_ROADMAP_20260929.md:
PHYS-001–007. Nie wykonuj PHYS-008; to checkpoint architekta.

Najpierw przeczytaj AGENTS.md i wskazane intake'y. Zachowaj aktualne
granice dowodów 327/337. Nie otwieraj holdoutu 335, nie zmieniaj
theta=2 ani Tprep=4 w protokole 333 i nie uruchamiaj dużych kampanii.

Pracuj w fin_physical_bridge_campaign_01/. Nie nadpisuj FIN son.
Najpierw N=1,2,3; N=4 tylko dla konkretnej potrzeby w kryterium.
Nie uruchamiaj domyślnego seed_checks.py. Najpierw napraw w kopii
kontrolę degeneracji/modów oraz fail-closed tests.

Celem nie jest nowa teoria wszystkiego. Celem jest:
dokładny słownik generalized Curie-Weiss,
A7-dependent statyczny fingerprint,
oddzielenie kinetyki od równowagi,
jawny source ledger
i wykonalny, kalibrowalny test wobec kontrmodelu.

Po każdym tasku oddaj REPORT.md, RESULTS.json, NONCONCLUSIONS.md,
REPLAY.md oraz hashe wejść. W RESULTS rozdziel execution_status
od scientific_verdict. Negatywny wynik zachowaj i zastosuj wskazaną
ścieżkę FAIL; nie dopisuj założenia po to, by otrzymać PASS.

Po PHYS-007 przygotuj wspólny HANDOFF.md, CLAIM_REGISTER.json,
NONCONCLUSIONS.md, MANIFEST.sha256 i minimalne polecenia replay.
Zatrzymaj się. Nie kupuj hardware, nie kontaktuj laboratorium,
nie pobieraj danych prywatnych i nie zmieniaj AGENTS.md bez
odrębnej zgody.

W HANDOFF odpowiedz: co jest standardową fizyką statystyczną,
co rzeczywiście zależy od FIN/A7, co ma niezależne uzasadnienie
fizyczne i co pozostało ręcznym założeniem.
~~~

### Stop conditions

1. Jeśli pozostaje tylko znana klasa mean-field z ręcznie wybraną A7,
   nazwij ją tak. Nie nazywaj tego fundamentalną teorią.
2. Jeśli realizacja tylko odtwarza zaprogramowane współczynniki,
   wynik dotyczy realizacji, nie pochodzenia A7.
3. Jeśli kontrast znika po dopuszczalnej kalibracji, odrzuć sondę lub
   zakres modelu; nie wprowadzaj nowego fitowanego parametru na każdy test.
4. Jeśli wybrana platforma nie realizuje rozkładu lub generatora,
   nie zastępuj ich podobieństwem funkcji celu.
5. Jeśli source lane potrzebuje nowej przesłanki, wróć do [S].
   Nie otwieraj ponownie starych no-go bez nowego typed atom.
6. Jeśli nie ma niezależnej predykcji i wszystkie pozytywne efekty są
   wspólne dla całej klasy, zatrzymaj inwestowanie w tę ścieżkę jako
   uzasadnienie fundamentalności. Dalsza praca aplikacyjna wymaga
   osobnego, jawnego celu.

## 11. Przegląd FIN son — stan na 2026-09-29

### Zakres

Zbadano całe dwa pliki: fin_core.py (234 linie) i seed_checks.py (85 linii).
Nie uruchomiono pełnego seed_checks. Wykonano ograniczone kontrole
N=2–3 oraz enumerację statyczną przy g bliskim zeru. Nie zmieniono plików autora.

### Co potwierdziły małe kontrole

- A7 ma rank 7, stałą przekątną, centrowanie i symetrię z residuals około 10^-15.
- Dla N=2 heat-bath, Metropolis i Barker zachowują podaną pi;
  błędy stationarity/detailed balance są około 10^-17 w tej arytmetyce.
- N=3 heat-bath daje rho=0.13143978619564534, zgodne z przyjętym fixture.
- Nowa kontrola słabego g odtwarza dokładny skończone-N wzór nachyleń S_k.
- To kontrole zgodności numerycznej, nie globalna certyfikacja ani
  dowód niezależności implementacji od każdego wcześniejszego wyboru modelu.

### Wykryte problemy

1. **slow_modes nie obsługuje ogólnej degeneracji.**
   W kontroli N=2,g=0 wystąpił KeyError: 4.
   Kod przypisuje całemu klastrowi eigenvalues tylko pierwszy k z
   rotacji. Multiplet może zawierać wiele irreps. Trzeba wyeksportować
   wszystkie sektory i ich residuals, nie ev[0].

2. **Symetryzacja nie może zastępować kontroli odwracalności.**
   Ss=(Ss+Ss.T)/2 wykonywane jest bez wymagania, by wejściowy defekt
   detailed balance był odpowiednio mały. Nieprawidłowy generator
   zostałby zastąpiony innym problemem spektralnym.

3. **V_d4_saddle jest lokalną procedurą, nie globalnym separatorem.**
   Warunek p_j=p_(4-j) nie nakłada nierówności P0=P1>=P2.
   Reflection-fixed subspace zawiera także localized states przy
   etykietach 2/8; sam constrained minimizer nie gwarantuje wyboru
   właściwego d4 saddle. Sprawdzaj root, indeks i właściwą half-wall.
   Nie używaj tej funkcji do zamknięcia Gamma=B4.

4. **g_eq używa heuristic branch discovery.**
   Brak znalezionego minimum zastępowany jest stałą 0.01, więc funkcja
   może być nieciągła. brentq nie certyfikuje globalnej coexistence.
   Wymagane są stationarity, końcowy energy residual i zakres gałęzi.

5. **Brak twardych gates optimizer success i saddle inertia.**
   Wynik BFGS/L-BFGS może być przyjęty mimo niepowodzenia. Parametry
   softmax mają też kierunek gauge, a V_d4 używa niestabilizowanego exp.

6. **Kontrmodel projektorowy ma inną normalizację.**
   Średnia czterech eigenvalues wynosi około 2.2004410061;
   trace-preserving średnia z multiplicities wynosi około 2.1801922868.
   Nie jest to automatycznie błąd, jeśli g/geq jawnie usuwa tę skalę,
   ale porównania przy tym samym g i przy tym samym g/geq nie są zamienne.

7. **Zakres „general A” wymaga ograniczenia.**
   Aktualne generatory są poprawne dla symmetric constant-diagonal A.
   Przy rozszerzeniu do niestałej przekątnej dodaj conditional diagonal
   correction albo jawnie zabroń tej klasy.

8. **G_FROZEN ma udokumentowane pochodzenie operacyjne.**
   R09 i wcześniejszy checker wiążą 5.145228719... z balansem
   k6/escape barriers. Komentarz „origin not documented” jest nieaktualny.
   Nie czyni to g stałą fizyczną.

9. **Domyślny runner jest szerszy niż jego nagłówek.**
   Mimo zapowiedzi N<=5 wykonuje static_S także dla N=6,8 i wiele
   optymalizacji g_eq. Zapis seed_results.json idzie do cwd, nie
   koniecznie obok skryptu. Nie wykonywać go jako bezkosztowego smoke testu.

### Werdykt dla FIN son

Użyteczna, czytelna implementacja zalążkowa i sensowny pomysł na rozdzielenie
static/dynamic fingerprints. Nie jest jeszcze bezpiecznym automatycznym
arbitrem alternatywnych modeli ani pakietem dowodowym. Najpierw PHYS-002,
potem wykorzystanie w nowej kampanii.

## 12. Zwięzłe wskazówki dla pozostawionych gałęzi 327/337

Nie są częścią pierwszej kampanii PHYS ani dodatkowo uruchomionymi taskami.

- 327: napraw dokładne enclosures wejść, LP i scalar minimum przed pełnym
  replay. Dla LP z sum p=1 i L_i<=p_i<=U_i dowolny próg z daje poprawne
  granice z + Σ min/max[(c_i-z)L_i,(c_i-z)U_i]. Dualny świadek nie wymaga
  uznania za dokładny zaokrąglonego greedy allocation. Wszystkie operacje
  nadal wymagają outward bounds. Dla dodatniej krzywizny obliczaj
  analityczne minimum kwadratu, nie wartość w zaokrąglonym argmin.
- 337: przy seedzie i polu na 0 działającą symetrią jednego protokołu jest
  stabilizator 0 (odbicie), nie dowolna translacja D12. Na D<=2 daje to
  1+6+36=43 orbit typów wejścia. Pełne D12 działa kowariantnie na rodzinie
  różnych seedów; nie utożsamiać ich w jednym kontrolowanym doświadczeniu.
- Wagi w obecnym low_defect_theta2_summary pochodzą z biased equilibrium
  conditioned on J=0. Pełny protokół 333 daje delta_seed exp(4 L_prep).
  Użyj rzeczywistych finite-time wag albo dolicz controller error.
- Retrospektywna diagnostyka bieżących tabel wskazuje, że weighted
  one-defect approximation ma TV error około 6.49e-5 dla N=7 i 1.66e-5
  dla N=10 na equilibrium-reference D<=2. Cancellation-free weighted
  pair bounds są odpowiednio około 6.75e-5 i 1.66e-5.
  Nie oznacza to zniknięcia pair effects ani transferu na nowe g.
- Ustalony D<=2 nie jest N-uniform przy stałym skończonym theta,g,Tprep>0.
  Każda kopia ma dodatnią, N-uniform szansę końcowego defektu; liczba
  defektów nie pozostaje ograniczona przy N→infinity. Wniosek o 78
  stanach jest finite-regime result.
- Lokalny transfer po g można wyprowadzać przez pochodne backward
  semigroup: u'=L_g u, v'=L_g v+(partial_g L_g)u. Nie otwieraj 335
  przed zamrożeniem radius/remainder i uwzględnieniem g-dependent
  basin/readout classification. Metodologiczny punkt odniesienia:
  https://eprints.maths.manchester.ac.uk/1218/
- Supplied NPZ są dostępne w bieżącym workspace, lecz są ignorowane
  przez Git. Handoff wykonawczy musi je dołączyć z hashami albo dać
  rzeczywisty producer; sam clone repo nie wystarczy.

## 13. Końcowa decyzja architektoniczna

Pierwszy most do znanej fizyki ma iść przez macierzowy model Curie–Weissa,
jego statyczny spectral fingerprint i kontrolowany sampler. Równolegle
należy uczciwie zakończyć source ledger A7. Nie ma obecnie podstaw,
by nazywać A7 fizycznie konieczną macierzą albo program FIN fundamentem QM/GR.

Jeżeli pierwsza kampania wykaże jedynie poprawną realizację standardowego
modelu z zadanymi couplings, należy to uznać za jej rzeczywisty wynik.
Dalsza inwestycja w fundamentalną interpretację wymaga nowego,
niezależnego źródła macierzy i przewidywania spoza kalibracji.

