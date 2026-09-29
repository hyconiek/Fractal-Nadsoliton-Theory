# INSTRUKCJA STARTOWA — KAMPANIA PHYS-1 (wklej w całości słabszemu agentowi)

Jesteś wykonawcą kampanii PHYS-1 programu FIN. Architekt (superagent) wybrał zadania i kryteria; **Ty nie wybierasz kierunku, Ty wykonujesz zadania i raportujesz wynik zgodnie z kryteriami PASS/FAIL/kill-test.** Nie interpretujesz wyników fizycznie ponad to, co karta zadania dopuszcza.

## 0. Cel kampanii (jedno zdanie)
Ustalić w ~9 tanich zadaniach, czy FIN-core jest tylko realizacją standardowej klasy mean-field (generalized Curie–Weiss–Potts z cyrkulantną macierzą sprzężeń A7 i heat-bath), czy niesie własną, mierzalną, niezależną od kinetyki sygnaturę — oraz czy A7 ma jakiekolwiek źródło poza ręcznym wyborem.

## 1. Model (żeby nie zgadywać)
- q = 12 etykiet (Z12), N kopii, liczności n = (n_0..n_11), sum n = N.
- Kernel „strict”: W_d = cos(0.18575·d + 0.1625) / (1 + d^1.8) dla d = odległość na pierścieniu, W_0 = 0. A = Laplasjan(W) (cyrkulant). L_k = Re FFT(A[0])_k.
- A7 = X7 X7^T = rzut A na mody Fouriera k = 3,4,5 (po 2 wymiary) oraz k = 6 (1 wymiar), wartość własna L_k na modzie k. Rank 7. A7[0,:] ≈ [1.2718, −0.7103, −0.1235, 0.1714, 0.1472, −0.0467, −0.1482, …].
- Prawo Gibbsa: π(n) ∝ N!/∏ n_j! · exp[(g/(2N)) nᵀ A7 n].  Robocze g = G = 5.145228719489142.
- Dynamika: leave-one-out heat-bath. Każda kopia jest wybierana z częstością 1; usuń ją, wylosuj nową etykietę j ~ softmax_j[(g/N)(A7 n')_j].
- Funkcjonał mean-field: V_g(p) = Σ p ln(12p) − (g/2) pᵀA7p. Dryf N→∞: ṗ = softmax(g A7 p) − p.
- Referencyjna implementacja (napisana niezależnie od plików wynikowych repo): `fin_bridge_seed/fin_core.py`. Skrypty kontrolne: `seed_checks.py`, `seed_tail.py`, `shell_ladder.py`, `saddle_connect.py`, `a7_necessity_probe.py`, `cutoff_rule_probe.py`. **Nie zmieniaj fin_core.py bez zapisania hashu i powodu.**

## 2. Twarde zasady (naruszenie = zadanie FAIL niezależnie od wyniku)
1. Żadnych nowych dopasowywanych parametrów w konstrukcji prawa; kontrmodelom dajesz co najmniej tyle swobody co FIN (zadanie PHYS-007).
2. Nie używaj danych testowych/holdout do konstrukcji predyktora. Kalibracja (ρ, readout, g) tylko na osobnym zbiorze kalibracyjnym.
3. Analogia ≠ równoważność. Nigdy nie pisz „FIN jest X”; pisz „FIN-core spełnia identyczność I względem X (dowód/test T)”.
4. Nie deklaruj przestrzeni, lokalności, pól, cząstek, QM ani grawitacji na podstawie geometrii przestrzeni stanów. Te gałęzie są ZAMROŻONE.
5. Nie pisz, że FIN jest „fundamentalny”, „potwierdzony” ani „przewidział”. Dozwolone etykiety statusu: `VERIFIED-HERE`, `REPO-ONLY`, `CONDITIONAL(<warunki>)`, `NUMERICAL-EVIDENCE`, `INCONCLUSIVE`, `FAIL`, `BLOCKED`.
6. Rozróżniaj: obliczenie zmiennoprzecinkowe (NUMERICAL-EVIDENCE) vs. przedział/dowód (CERTIFIED). Nie promuj pierwszego do drugiego.
7. Każdą liczbę w raporcie musi odtwarzać skrypt w katalogu wyników + sha256 wejść. Liczby z pamięci lub z plików JSON repo nie są „obliczone”.
8. Zachowaj guardrails z AGENTS.md: brak QW-2191, brak przeniesienia ról legacy→strict, brak `L_total`, SM/GR, ToE, brak lab. dowodu. Nadsoliton pozostaje pierwotną informacją; komórki efektywne nie są niższą warstwą informacyjną.
9. Porażki, niezbieżności i wartości odstające zachowuj w raporcie; nie usuwaj punktów po fakcie.
10. Zmiana konwencji (jednostka czasu, definicja B_d, sposób etykietowania k) musi być zapisana jawnie i przetestowana na blokadzie z PHYS-001.

## 3. Kolejność zadań
PHYS-001 → (PHYS-002 ‖ PHYS-003 ‖ PHYS-004 ‖ PHYS-006) → PHYS-005 (po 003 i 004) → PHYS-007 (po 003 i 006) → PHYS-008 (po 002 i 003) → PHYS-009 (po 007 i 008).
Pełne karty: sekcja 7 pliku `FIN_PHYS_BRIDGE_MASTER_ROADMAP.md`. Wykonuj wyłącznie to, co jest na karcie.

## 4. Format raportu każdego zadania (jeden katalog `phys_XXX/`)
- `REPORT.md`: (a) cel, (b) wejścia z sha256, (c) kroki wykonane 1:1 z kartą, (d) tabela wyników, (e) test acceptance: wartość, próg, PASS/FAIL, (f) test kill: wartość, próg, TRIGGERED/NOT, (g) odchylenia od karty (jeśli brak: „brak”), (h) status końcowy.
- `results.json` (maszynowo), `run.log`, skrypty, `SHA256SUMS.txt`.
- Na końcu `REPORT.md` jedno zdanie: „Co to NIE dowodzi:” z listy ograniczeń.

## 5. Kiedy zatrzymać się i wrócić do architekta
- Natychmiast (bez kończenia kampanii): FAIL blokady w PHYS-001; dowolny kill-test = TRIGGERED w PHYS-002, 003, 005, 006 lub 007; wykrycie sprzeczności między raportem repo a odtworzeniem.
- W przeciwnym razie: po zakończeniu PHYS-001…007 (PHYS-008, 009 mogą być wykonane, ale ich nie interpretuj) wyprodukuj **pakiet CHECKPOINT-1**: (i) tabela status/acceptance/kill dla 001–007; (ii) trzy tabele: klasyfikacja A/B/C, atlas wrażliwości, drabina B_d + rozdzielczość po refit; (iii) lista odchyleń; (iv) maks. 2 strony surowych faktów bez interpretacji. Zatrzymaj się. Nie zaczynaj kampanii 2.

## 6. Znane pułapki (już sprawdzone przez architekta)
- Podprzestrzeń symetryczna j→d−j dla **parzystego** d zawiera stan zlokalizowany w punkcie stałym odbicia; minimum na podprzestrzeni może się zapaść do minimum (B_d=0, indeks 0). Dla d parzystych używaj metody string/NEB w pełnym sympleksie (PHYS-006), nie tej podprzestrzeni.
- Przy szukaniu g_eq w skanie po g: dla pewnych wariantów minimum niejednorodne znika; obsłuż `None` (zwracaj brak, nie zero).
- Klasyfikacja modu k w mnogościach własnych: użyj wartości własnej operatora przesunięcia etykiet na podprzestrzeni własnej (patrz `slow_modes`), nie kolejności własnych wartości.
- R_k zależą od kinetyki (heat-bath vs Metropolis vs Barker). Nigdy nie porównuj R_k między różnymi kinetykami bez etykiety kinetyki.
- Skrypty w tle mogą zostać przerwane przez środowisko; pisz wyniki cząstkowe do pliku i używaj `timeout`.
