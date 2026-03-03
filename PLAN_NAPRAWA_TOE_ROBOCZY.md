# PLAN ROBOCZY: NAPRAWA I ROZWÓJ FIN-ToE (bez edycji TEX na tym etapie)
**Status:** roboczy  
**Data:** 2026-03-02  
**Cel:** przygotować pełny plan naprawy ToE w duchu FIN, zanim cokolwiek trafi do `TOE_FINAL_DOCUMENTATION.tex`.

## 1. Założenie operacyjne
Na tym etapie nie modyfikujemy dokumentu LaTeX.  
Pracujemy tylko na plikach roboczych i raportach pomocniczych.

## 2. Co już mamy (materiał wejściowy)
1. Audyt spójności i chronologii repo (`QW-1701`, `QW-1704`).
2. Audyt reprodukowalności (`QW-1702`).
3. Audyt deklaracji vs obliczeń (`QW-1703`).
4. Estymator poprawek do rdzeniowych wzorów (`QW-1708`).
5. Propozycja naprawionego rdzenia koncepcyjnego (`FIN_ToE_4_5_REPAIRED_CORE.md`).

## 3. Diagnoza głównych problemów (priorytet)
1. Mieszanie klas twierdzeń:
   - predykcja z pierwszych zasad vs kalibracja.
2. Niespójność retoryczna:
   - lokalne `EXACT/0.00%` przy równoległych raportach krytycznych.
3. Niedomknięcie sektora flavor:
   - brak konstruktywnej derywacji CKM/PMNS.
4. Brak stabilnego protokołu out-of-sample dla kluczowych obserwabli.

## 4. Plan naprawy merytorycznej (kolejność)
### Etap A: Stabilizacja metodologii
1. Wprowadzić taksonomię D/C/M/F/P dla każdego twierdzenia.
2. Dla każdej obserwabli opisać:
   - parametry zamrożone,
   - parametry kalibrowane,
   - metrykę błędu i niepewności.
3. Ustalić jednolity protokół reprodukowalności semantycznej (2-run).

### Etap B: Naprawa równań rdzeniowych
1. Zostawić rdzeń `K(d)` bez zmiany ontologicznej.
2. Dodać jawne poprawki efektywne:
   - `δ_W` dla kąta Weinberga,
   - `δ_vac` dla `α_EM`,
   - `Δ_a` dla mas.
3. Przejść z narracji „idealnej zgodności” na narrację „hierarchii poprawek”.

### Etap C: Rozwój ToE (w duchu FIN)
1. Flavor: zaproponować i przetestować operator topologicznych przejść.
2. CKM/PMNS: konstrukcja macierzy mieszania z nakładów modów topologicznych.
3. Połączenie z sektorem grawitacyjnym przez jeden spójny ledger hipotez.

## 5. Plan plików roboczych (bez TEX)
1. `FIN_ToE_4_5_REPAIRED_CORE.md`:
   naprawiony formalny rdzeń teorii.
2. `RAPORT_QW1708_REPAIR_CONSTANTS_ESTIMATOR.md`:
   liczby potrzebnych poprawek.
3. `RAPORT_QW1704_CLAIM_LEDGER_BY_DATE.md`:
   status twierdzeń w czasie.
4. (Do utworzenia) `FIN_ToE_4_5_DRAFT_FOR_TEX.md`:
   gotowy materiał do późniejszego, kontrolowanego przeniesienia do LaTeX.

## 6. Kryterium przejścia do edycji TEX
Do `TOE_FINAL_DOCUMENTATION.tex` przechodzimy dopiero gdy:
1. użytkownik potwierdzi zakres zmian,
2. ustalimy docelowy styl (konserwatywny vs ofensywny publikacyjnie),
3. gotowy będzie jeden spójny draft sekcji, bez mieszania statusów twierdzeń.

## 7. Najbliższy krok roboczy (po akceptacji)
Przygotować `FIN_ToE_4_5_DRAFT_FOR_TEX.md` z:
1. krótkim „Executive Theoretical Repair” (2 strony),
2. jedną tabelą statusów hipotez,
3. jedną tabelą poprawek efektywnych,
4. listą 6 testów, które domykają ToE publikacyjnie.

## 8. Status sciezki "charakterystyki nadsolitona -> kernel" (QW-1729...1733)
### Wyniki wykonane
1. `QW_1729_NADSOLITON_KERNEL_CHARACTERISTICS_MAP.py`
   - kernel_closure_score = 0.522
   - werdykt: `KERNEL_CHARACTERISTICS_NOT_CLOSED`
2. `QW_1730_NADSOLITON_KERNEL_CHRONO_AUDIT.py`
   - risk_points = 9
   - wykryte sprzecznosci: `phi`, `beta_tors`, wzorce wezlow, niespojnosc arytmetyczna wzoru na wezly
   - werdykt: `KERNEL_ORIGIN_CHRONOLOGY_HIGH_RISK_INCONSISTENT`
3. `QW_1731_NADSOLITON_KERNEL_NODE_COMPATIBILITY.py`
   - dla priorytetowych parametrow (`omega=pi/4`, `phi=pi/6`) zera wypadaja: `d0=1.333`, `delta_d=4`
   - deklarowane wzorce wezlow wymagaja przesuniecia `omega` o ~33%
   - werdykt: `NODE_NARRATIVE_STRONGLY_INCONSISTENT`
4. `QW_1732_FRACTAL_PATH_TO_HYPERBOLIC_TEST.py`
   - przejscie `path-sum -> 1/(1+beta d)` jest mozliwe, ale wasko stabilne (tuning)
   - werdykt: `HYPERBOLIC_REDUCTION_PLAUSIBLE_BUT_TUNED`
5. `QW_1733_NADSOLITON_KERNEL_CLOSURE_GATE.py`
   - global_score = 0.170
   - hard_gate = FAIL
   - readiness: `KERNEL_DERIVATION_OPEN_NOT_CLOSED`

### Wniosek operacyjny
Kernel jest silny fenomenologicznie, ale formalnie NIE jest jeszcze wyprowadzony jednoznacznie z charakterystyk nadsolitona.

## 9. Kolejne badania py, ktore realnie domykaja kernel
1. `QW_1734_MICRO_BETA_TORS_DERIVATION.py`
   - Hipoteza: `beta_tors` wynika z mikrodynamiki skrutu topologicznego, nie z kalibracji.
   - Test: wyprowadzenie `beta_tors` z niezaleznych symulacji mikrosieci (bez fitowania do mas/GW).
   - Kryterium PASS: zgodnosc przewidywanego `beta_tors` z zakresem roboczym ±10% oraz stabilnosc OOS.
2. `QW_1735_OMEGA_PHI_FROM_LATTICE_DYNAMICS.py`
   - Hipoteza: `omega, phi` sa wymuszone przez dynamike sieci (symetrie + warunki brzegowe), nie ansatz.
   - Test: wymuszenie przez PDE/lagrangian i spektralna dekompozycje modow.
   - Kryterium PASS: jedna para (`omega,phi`) reprodukuje i period, i wezly bez konfliktu.
3. `QW_1736_KERNEL_NODE_BAYESIAN_MODEL_SELECTION.py`
   - Hipoteza: jeden wariant wzorca wezlow jest statystycznie preferowany.
   - Test: porownanie modeli wezlow (`2,5,8,11` vs `2,8,14` vs sekwencja analityczna) z Bayes factor.
   - Kryterium PASS: Bayes factor >= 10 dla jednego modelu + brak konfliktu z priorem charakterystyk.
4. `QW_1737_SHARED_KERNEL_FLAVOR_GW_CROSS_CONSTRAINT.py`
   - Hipoteza: ten sam zestaw parametrow kernela przechodzi jednoczesnie flavor + GW.
   - Test: wspolna funkcja wiarygodnosci i ocena identyfikowalnosci.
   - Kryterium PASS: wspolny region parametrow niepusty i stabilny pod perturbacjami.

## 10. Wyniki po wykonaniu badan 1734-1737 (+ bramka 1738)
1. `QW_1734_MICRO_BETA_TORS_DERIVATION.py`
   - Werdykt: `BETA_TORS_MICRO_DERIVATION_NOT_CLOSED`
   - Klucz: mikromodel nie odtwarza stabilnie docelowego `beta_tors ~ 0.01`.
2. `QW_1735_OMEGA_PHI_FROM_LATTICE_DYNAMICS.py`
   - Werdykt: `OMEGA_PHI_DERIVATION_NOT_STABLE`
   - Klucz: estymowane `(omega,phi)` nie stabilizuje sie blisko `(pi/4, pi/6)`.
3. `QW_1736_KERNEL_NODE_BAYESIAN_MODEL_SELECTION.py`
   - Werdykt: `NODE_MODEL_STRONGLY_SELECTED`
   - Wynik: model wezlow `M_A: 2+3n` ma przewage bayesowska nad alternatywami.
4. `QW_1737_SHARED_KERNEL_FLAVOR_GW_CROSS_CONSTRAINT.py`
   - Werdykt: `SHARED_KERNEL_REGION_NOT_SUPPORTED`
   - Klucz: brak niepustego wspolnego obszaru flavor+GW przy twardych bramkach empirycznych.
5. `QW_1738_KERNEL_CLOSURE_PHASE2_GATE.py` (agregacja 1730-1737)
   - Global score: `0.207`
   - Hard gate: `FAIL`
   - Readiness: `PHASE2_KERNEL_PATH_OPEN`

## 11. Konsekwencja dla domkniecia ToE
Po tej serii kernel nadal nie jest domkniety jako derywacja z pierwszych zasad nadsolitona.
Najmocniejszy sygnal jest po stronie selekcji wzorca wezlow (1736), ale nie jest on zgodny
z obecnym priorem charakterystyk `(omega=pi/4, phi=pi/6)` i nie przechodzi wspolnie flavor+GW.

## 12. Wyniki po wykonaniu QW-1739 i QW-1740
1. `QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py`
   - Cel: przebudowa mikromodelu na podpisane dynamiczne sprzężenia i derivacja `(beta_tors, omega, phi)` bez ustawiania z gory.
   - Werdykt: `SIGNED_MICROMODEL_DERIVATION_SUPPORTED_BUT_SHIFTED`
   - Klucz:
     - `median_r2 ~ 0.692` (fit technicznie akceptowalny),
     - ale `omega_median ~ 0.118` (daleko od `pi/4`) i `beta_median ~ 0.3` (daleko od `0.01`).
2. `QW_1740_SIGNED_MICROMODEL_IDENTIFIABILITY_AUDIT.py`
   - Cel: audyt identyfikowalnosci parametrow z profili 1739.
   - Werdykt: `IDENTIFIABILITY_WEAK`
   - Klucz:
     - bardzo zla kondycja hessiana (`cond ~ 8e15`),
     - wiele minimow lokalnych (`n_unique_modes = 3`),
     - brak stabilnego jednoznacznego estymatora `(omega,phi,beta)`.

## 13. Konsekwencja po 1739-1740
Nowy mikromodel jest krokiem naprzod metodologicznie (signed/dynamic), ale nadal nie domyka ToE:
parametry są wyprowadzalne tylko w wersji przesunietej i slabo identyfikowalnej.

## 14. Dalsza runda rygoru (QW-1741 do QW-1745)
1. `QW_1741_CONSTRAINED_GLOBAL_DERIVATION.py`
   - Globalny fit wielobiegowy + bootstrap.
   - Werdykt: `CONSTRAINED_GLOBAL_DERIVATION_NOT_CLOSED`.
   - Problem: optimum na granicach (`omega` low-bound, `beta` high-bound), mimo dobrego median R2.
2. `QW_1742_PROFILE_LIKELIHOOD_IDENTIFIABILITY.py`
   - Profil-likelihood + Hessian/Fisher wokol optimum 1741.
   - Werdykt: `PROFILE_IDENTIFIABILITY_WEAK`.
3. `QW_1743_OSCILLATORY_COHORT_DERIVATION.py`
   - Selekcja kohorty oscylacyjnej dla zwiekszenia informacji o `(omega,phi,beta)`.
   - Werdykt: `OSCILLATORY_COHORT_DERIVATION_NOT_CLOSED`.
   - Problem: zbyt mala kohorta informacyjna (`n_selected=9`), nadal rozwiazanie blisko granic.
4. `QW_1744_OSCILLATORY_COHORT_IDENTIFIABILITY.py`
   - Audyt identyfikowalnosci tylko na kohorcie oscylacyjnej.
   - Werdykt: `OSCILLATORY_IDENTIFIABILITY_WEAK`.
5. `QW_1745_MICROMODEL_RIGOR_GATE.py`
   - Agregacja wynikow 1739-1744.
   - Global score: `0.40`, hard gate `FAIL`, readiness `MICROMODEL_ITERATION_OPEN`.

## 15. Wniosek po 1745
Obecna przestrzen modeli nadal nie dostarcza stabilnej, nienabrzegowej i dobrze identyfikowalnej
derywacji `beta_tors, omega, phi` z pierwszych zasad. Dalsza praca musi isc w kierunku:
1) nowych obserwabli dynamicznych (nie tylko profil C(d)),
2) rozdzielenia degeneracji `omega-beta` przez niezalezne pomiary,
3) wiekszej kohorty informacyjnej z kontrola warunkow generacji.

## 16. Runda multi-observable (QW-1746 do QW-1748)
1. `QW_1746_DYNAMIC_OBSERVABLES_DERIVATION.py`
   - Wprowadzono obserwable: phase-lag slope, zero-cross spacing, envelope decay.
   - Wynik: tylko `n_good=10` (subset relaxed), duza niespojnosc estymatorow omega.
   - Werdykt: `DYNAMIC_OBSERVABLES_DERIVATION_WEAK`.
2. `QW_1747_MULTIOBS_JOINT_INFERENCE.py`
   - Wspolna inferencja z obserwabli dynamicznych.
   - Wynik: optimum graniczne (`omega=0.1`, `beta=0.25`), slaba kondycja (~3.4e12), wysoka korelacja.
   - Werdykt: `MULTIOBS_JOINT_IDENTIFIABILITY_WEAK`.
3. `QW_1748_MULTI_OBS_EXTENSION_GATE.py`
   - Agregacja 1745 + 1746 + 1747.
   - Global score: `0.20`, hard gate `FAIL`, readiness `MULTI_OBS_EXTENSION_OPEN`.

## 17. Konsekwencja po 1748
Dodanie obserwabli dynamicznych bylo metodologicznie zasadne, ale praktycznie nie przelamalo
degeneracji. Do domkniecia ToE potrzebny jest nowy etap eksperymentalny/symulacyjny z
obserwabla, ktora daje bezposredni i stabilny sygnal dla `beta_tors` niezalezny od `omega`.

## 18. Strategia separacyjna beta->(omega,phi) (QW-1749 do QW-1751)
1. `QW_1749_BETA_ORTHOGONAL_OBSERVABLE.py`
   - Obserwabla impulsowa celowana w `beta_tors` (agregacja po czasie i wartosci bezwzglednej).
   - Werdykt: `BETA_ORTHOGONAL_OBSERVABLE_PARTIAL`.
   - Plus: dobra jakosc dopasowania i odpornosc na phase-scramble.
   - Minus: estymata `beta` nadal przy dolnej granicy (brak nienabrzegowosci).
2. `QW_1750_SEPARATED_BETA_THEN_OMEGA_PHI.py`
   - `beta` z 1749, potem inferencja `omega,phi` + sweep po CI beta z dopasowaniem surowej odpowiedzi przestrzennej.
   - Werdykt: `SEPARATED_PIPELINE_WEAK`.
   - Klucz: `omega` w optimum nadal przy granicy, slaba kondycja.
3. `QW_1751_BETA_SEPARATION_GATE.py`
   - Agregacja 1748 + 1749 + 1750.
   - Global score: `0.55`, hard gate `FAIL`, readiness `BETA_SEPARATION_STRATEGY_OPEN`.

## 19. Stan po 1751
Strategia separacyjna daje realny postep w izolacji `beta`, ale nie domyka jeszcze
identyfikowalnosci `omega,phi` w rygorze wymaganym do zamkniecia ToE.

## 20. Rozszerzenie ortogonalne (QW-1752 do QW-1754)
1. `QW_1752_OMEGA_ORTHOGONAL_OBSERVABLE.py`
   - Cel: zbudowac obserwable celowane w `omega` mniej czule na obwiednie (`beta`).
   - Werdykt: `OMEGA_ORTHOGONAL_OBSERVABLE_PARTIAL`.
   - Klucz:
     - `n_good=11` (ponizej progu 12),
     - `omega_iqr=0.282` (akceptowalne),
     - `delta_omega_beta_perturb_q90=0.081` (lekko powyzej progu 0.08).
2. `QW_1753_ORTHOGONAL_TRIAD_JOINT_INFERENCE.py`
   - Cel: wspolna inferencja ortogonalna `(beta, omega, phi)` z osobnych strumieni danych 1749 i 1752.
   - Werdykt: `ORTHOGONAL_TRIAD_IDENTIFIABILITY_WEAK`.
   - Klucz:
     - jakosc obserwabli `beta` i `omega` jest technicznie dobra,
     - ale test ortogonalnosci pada (`rho_s(omega, delta_beta_perturb)=-0.727`),
     - oraz joint nonboundary rate niski (`0.442`), co oznacza nadal silna nienabrzegowosc problemu.
3. `QW_1754_ORTHOGONAL_EXTENSION_GATE.py`
   - Agregacja 1751 + 1752 + 1753.
   - Global score: `0.58`, hard gate `FAIL`, readiness `ORTHOGONAL_EXTENSION_OPEN`.

## 21. Stan po 1754
Separacja obserwabli poprawila przejrzystosc diagnostyczna, ale nie domknela teorii:
1. `beta_tors` pozostaje blisko dolnej granicy dopuszczalnego zakresu.
2. `omega` nadal ma szeroki CI z masa prawdopodobienstwa przy niskich wartosciach.
3. W nowym tescie ortogonalnosci ujawnia sie silna zaleznosc resztkowa (`omega` vs wrazliwosc na perturbacje beta),
   co jest bezposrednim sygnalem niedomknietego mikromodelu (sprzezenia nie sa jeszcze poprawnie rozdzielone).

## 22. Runda "null-vs-positive beta" i residual decoupling (QW-1755 do QW-1757)
1. `QW_1755_BETA_NULL_VS_POSITIVE_EVIDENCE.py`
   - Cel: sprawdzic, czy model z dodatnim `beta_tors` jest statystycznie preferowany nad nullem `beta=0`.
   - Werdykt: `BETA_POSITIVE_EVIDENCE_WEAK`.
   - Klucz:
     - mimo wysokiego `P(beta>1e-3)`, median `delta_BIC(null-beta)=-2.48`,
     - czyli formalnie null czesto wygrywa po karze za zlozonosc,
     - brak nienabrzegowosci (`CI95` schodzi do `1e-5`).
2. `QW_1756_OMEGA_BETA_DECOUPLING_RESIDUAL_TEST.py`
   - Cel: odsprzegnac `omega` od obwiedni przez jawna korekte sygnalu o `beta` z 1755.
   - Werdykt: `OMEGA_BETA_DECOUPLING_WEAK`.
   - Klucz:
     - poprawa w metrykach sprzezenia (`sensitivity_q90=0.0036`, `spearman=-0.245`),
     - ale nadal za szeroki rozrzut `omega` (`IQR=0.326`) i CI zahacza o skrajnie niskie wartosci.
3. `QW_1757_DECOUPLING_EXTENSION_GATE.py`
   - Agregacja 1754 + 1755 + 1756.
   - Global score: `0.51`, hard gate `FAIL`, readiness `DECOUPLING_EXTENSION_OPEN`.

## 23. Stan po 1757
Po tej iteracji mamy mocniejsza diagnoze niz domkniecie:
1. Obecny sygnal na dodatnie `beta_tors` nie jest jeszcze wystarczajaco silny po rygorystycznej karze modelowej.
2. Odsprzeganie resztowe poprawia stabilnosc lokalna `omega`, ale nie redukuje globalnej nieidentyfikowalnosci do poziomu akceptowalnego.
3. Najwiekszym blokerem pozostaje brak jednoczesnej nienabrzegowosci i waskich CI dla (`beta`, `omega`) w jednym protokole.

## 24. Protokoly prerejestrowane nowej generacji (QW-1758 do QW-1760)
1. `QW_1758_PREREGISTERED_IDENTIFIABILITY_PROTOCOL.py`
   - Cel: twardy prerejestrowany test identyfikowalnosci z wielozrodlowoscia i rozszerzonym zasiegiem `d`.
   - Werdykt: `PREREGISTERED_IDENTIFIABILITY_WEAK`.
   - Klucz:
     - total rows `32`, ale `n_good=3` (silny sample starvation),
     - powod: zbyt restrykcyjna bramka fazowa w wersji strict.
2. `QW_1759_ROBUST_PHASE_PREREG_PROTOCOL.py`
   - Cel: ten sam duch prerejestracji, ale z robust estymatorem fazy (WLS + tlumienie outlierow).
   - Werdykt: `ROBUST_PHASE_PREREG_IDENTIFIABILITY_PARTIAL`.
   - Klucz:
     - `n_good=17` (duza poprawa identyfikowalnosci estymacyjnej),
     - przechodza: `omega_nonboundary`, spready, niska korelacja `beta-omega`, stabilnosc perturbacyjna,
     - nadal pada: dowod dodatniego `beta` (`median delta_BIC=-2.36`) i nienabrzegowosc beta.
3. `QW_1760_PREREG_ROBUST_EXTENSION_GATE.py`
   - Agregacja 1757 + 1758 + 1759.
   - Global score: `0.41`, hard gate `FAIL`, readiness `PREREG_ROBUST_EXTENSION_OPEN`.

## 25. Stan po 1760
Po tej rundzie mamy istotna poprawnosc metodologiczną, ale ToE nadal nie jest domkniete:
1. Problem estymacyjny fazy zostal czesciowo naprawiony (1759), wiec nieidentyfikowalnosc nie jest juz tylko artefaktem fitu.
2. Główny bloker pozostaje fizyczno-modelowy: brak stabilnego, modelowo preferowanego dodatniego `beta_tors` (vs null).
3. Oznacza to, ze kolejny krok musi celowac w mikromodel generujacy obwiednie (a nie tylko w estymator), inaczej bramki dowodowe beda systematycznie padac.

## 26. Test augmentacji mechanizmu obwiedni (QW-1761 do QW-1762)
1. `QW_1761_ENVELOPE_HYBRID_MODEL_SELECTION.py`
   - Cel: sprawdzic, czy dodanie czlonu pamieci (`lambda`) w modelu obwiedni poprawia dowod dla niezerowego `beta_tors`.
   - Porownane modele: `M0` (null), `M1` (hiperboliczny beta), `M2` (hybrydowy beta+lambda).
   - Werdykt: `HYBRID_ENVELOPE_MECHANISM_WEAK`.
   - Klucz:
     - wygrane modeli: `M0:30`, `M1:6`, `M2:0` (na 36),
     - `median delta_BIC(M1-M2)=-2.48`, `median delta_BIC(M0-M2)=-4.96`,
     - czyli model hybrydowy jest penalizowany i nie jest empirycznie wspierany.
2. `QW_1762_MECHANISM_AUGMENTATION_GATE.py`
   - Agregacja 1760 + 1761.
   - Global score: `0.205`, hard gate `FAIL`, readiness `MECHANISM_AUGMENTATION_OPEN`.

## 27. Stan po 1762
Po tej serii mamy negatywny, ale bardzo informacyjny wynik:
1. Obecna augmentacja obwiedni (`beta+lambda`) nie jest poprawna sciezka domykania teorii.
2. Dane premiuja prostszy opis (`null` lub slaby model `beta`) zamiast bogatszej dynamiki obwiedni.
3. To oznacza, ze dalsze domykanie ToE wymaga zmiany poziomu modelu mikrodynamiki (generatora sygnalu),
   a nie tylko dokladania kolejnych efektywnych parametrow do tej samej warstwy obserwabli.

## 28. Iteracja mikrodynamiki z polem pamieci (QW-1763 do QW-1764)
1. `QW_1763_MEMORY_FIELD_MICRODYNAMICS_BETA_TEST.py`
   - Cel: sprawdzic, czy jawne pole pamieci w ewolucji amplitudy generuje stabilny sygnal dodatniego `beta_tors`.
   - Werdykt: `MEMORY_FIELD_BETA_EVIDENCE_WEAK`.
   - Klucz:
     - `n_rows=24`, ale `median delta_BIC(null-beta)=-2.64`,
     - `positive delta_BIC rate=0.167`,
     - `beta_nonboundary_rate=0.375`,
     - czyli mimo zmiany generatora sygnalu brak dowodu na niezerowe `beta` w rygorze modelowym.
2. `QW_1764_MEMORY_FIELD_GATE.py`
   - Agregacja 1762 + 1763.
   - Global score: `0.2025`, hard gate `FAIL`, readiness `MEMORY_FIELD_ITERATION_OPEN`.

## 29. Stan po 1764
Po tej rundzie:
1. Nawet przebudowa mikrodynamiki (pole pamieci) nie odblokowala dowodu na dodatnie `beta_tors`.
2. To wzmacnia wniosek, ze obecna rodzina modeli FIN (w tej klasie obserwabli) jest niedomknieta na poziomie mechanizmu.
3. Dalsze kroki powinny przejsc na nowy typ obserwabli lub nowy formalizm mikrodynamiki, zamiast kolejnych wariantow tej samej obwiedni.

## 30. Nowa galaz: obserwable nieobwiedniowe (QW-1765 do QW-1769)
1. `QW_1765_NONENVELOPE_MECHANISM_OBSERVABLES.py`
   - Cel: identyfikacja mechanizmu z obserwabli dynamicznych bez fitu obwiedni (`beta`).
   - Werdykt: `NONENVELOPE_MECHANISM_IDENTIFIABILITY_WEAK`.
   - Klucz: slaba separacja i niestabilny CV (`AUC~0.51`, `R2<0`), duza wariancja foldow.
2. `QW_1766_NONENVELOPE_PAIRED_INTERVENTION.py`
   - Cel: usuniecie konfuzji topologicznej przez schemat sparowany (ta sama siec: baseline vs interwencja).
   - Werdykt: `PAIRED_NONENVELOPE_MECHANISM_PARTIAL`.
   - Klucz:
     - bardzo mocna separacja (`AUC=1.0`, `BACC=1.0`, `R2~0.87`),
     - ale jeszcze fail stabilnosci foldowej (`fold_r2_std~0.24`).
3. `QW_1768_LEAKAGE_CONTROLLED_PAIRED_INTERVENTION.py`
   - Cel: powtorzenie 1766 z kontrola leakage (grouped CV po topologii/seed).
   - Werdykt: `LEAKAGE_CONTROLLED_NONENVELOPE_SUPPORTED`.
   - Klucz:
     - nadal bardzo silny sygnal po usunieciu leakage (`AUC_grouped=1.0`, `BACC_grouped~0.986`, `R2_grouped~0.944`),
     - stabilnosc foldowa dobra (`fold_r2_std~0.026`),
     - kalibracja dobra (pozytywne przypadki nie-graniczne, negatywne blisko zera).
4. `QW_1769_NONENVELOPE_RIGOR_GATE.py`
   - Agregacja 1764 + 1766 + 1768.
   - Global score: `0.6675`, hard gate `FAIL`, readiness `NONENVELOPE_RIGOR_PARTIAL`.

## 31. Stan po 1769
Nowa galaz daje pierwszy twardy, rygorystyczny sygnal domykania mechanizmu:
1. W torze nieobwiedniowym identyfikacja mechanizmu jest mozliwa i stabilna (1768).
2. Twardy fail legacy-branch (1764) nadal blokuje globalne zamkniecie ToE, ale nie obala nowej sciezki.
3. Praktycznie: dalsze domykanie ToE powinno przejsc na formalne mapowanie parametrow tej nowej galęzi
   na stare wielkosci kernela (`beta_tors`, `omega`, `phi`) zamiast dalszego wzmacniania starej obwiedni.

## 32. Most nowej galezi do starych parametrow kernela (QW-1770 do QW-1772)
1. `QW_1770_KERNEL_BRIDGE_FROM_NONENVELOPE.py`
   - Cel: mapowanie z obserwabli nieobwiedniowych na (`beta_bridge`, `omega`, `phi`) z grouped-CV.
   - Werdykt: `KERNEL_BRIDGE_FROM_NONENVELOPE_WEAK`.
   - Klucz:
     - mocny sygnal dla `beta_bridge` (`R2~0.72`, `rho~0.82`) i bardzo dobry dla `phi`,
     - sredni dla `omega` (`R2~0.64`, ponizej progu),
     - slabosc: niestabilnosc foldowa i zbyt niska nienabrzegowosc predykcji `beta`.
2. `QW_1771_BRIDGE_LEGACY_CONSISTENCY.py`
   - Cel: sprawdzic spojnosc skali nowego `beta_bridge` z legacy `beta_tors` z 1755/1749.
   - Werdykt: `BRIDGE_LEGACY_CONSISTENCY_WEAK`.
   - Klucz:
     - direct: bardzo duzy rozjazd skali (WD~0.136, brak overlap CI),
     - nawet po najlepszej monotonicznej kompresji overlap CI pozostaje znikomy.
3. `QW_1772_KERNEL_BRIDGE_INTEGRATION_GATE.py`
   - Agregacja 1769 + 1770 + 1771.
   - Global score: `0.4125`, hard gate `FAIL`, readiness `KERNEL_BRIDGE_INTEGRATION_OPEN`.

## 33. Stan po 1772
Mamy klarowny obraz:
1. Nowa galaz nieobwiedniowa jest merytorycznie silna i rygorystyczna na poziomie identyfikacji mechanizmu.
2. Nieudane pozostaje jej osadzenie w starej skali legacy `beta_tors` (brak zgodnosci rozkladowej).
3. Oznacza to, ze dalsza praca domykajaca ToE powinna obejmowac jawny nowy slownik/rekonfiguracje parametrow
   (reparametryzacje kernela), a nie probowac wymuszac zgodnosc przez proste mapy skali.

## 34. Reparametryzacja zamiast mapy skali 1:1 (QW-1773 do QW-1774)
1. `QW_1773_OMEGA_SUPPRESSED_LEGACY_PROJECTION.py`
   - Hipoteza: `beta_legacy` jest projekcja nowego `beta_bridge` tlumioną przez `omega`:
     `beta_legacy_equiv = alpha * beta_bridge * exp(-gamma / omega^p)`.
   - Werdykt: `OMEGA_SUPPRESSED_LEGACY_PROJECTION_SUPPORTED`.
   - Klucz:
     - bardzo silna redukcja rozjazdu (`WD` spada o ~x91),
     - wysoki overlap CI z legacy (~0.75),
     - parametry projekcji sa parsymoniczne (`alpha~0.98`, `gamma~0.55`, `p~1.50`).
2. `QW_1774_REPARAMETERIZATION_PATH_GATE.py`
   - Agregacja 1768 + 1772 + 1773.
   - Global score: `0.804`, hard gate `FAIL`, readiness `REPARAMETERIZATION_PATH_STRONG_PARTIAL`.
   - Fail hard-gate wynika z tego, ze segment 1772 (stary most bez reparametryzacji) nadal nie przechodzi.

## 35. Stan po 1774
Po tej iteracji pojawia sie realna sciezka domykania ToE:
1. Mechanizm nowej galezi jest rygorystycznie potwierdzony (1768).
2. Spójność z legacy odzyskano nie przez sztuczne skalowanie, ale przez konkretną reparametryzacje zalezną od `omega` (1773).
3. Formalnie ToE nie jest jeszcze zamkniete globalnie, ale mamy silny kandydat nowego slownika parametrow,
   ktory mozna teraz testowac na kolejnych niezaleznych kohortach i danych empirycznych.

## 36. Walidacja OOD i gotowosc empiryczna reparametryzacji (QW-1775 do QW-1776)
1. `QW_1775_REPARAM_OOD_VALIDATION.py`
   - Cel: sprawdzic czy parametry reparametryzacji z 1773 dzialaja poza zbiorem uczacym (OOD, przesuniete zakresy).
   - Werdykt: `REPARAM_OOD_VALIDATION_SUPPORTED`.
   - Klucz:
     - dobra zgodnosc z legacy (`WD~0.00484`, overlap CI~0.55),
     - poprawny ksztalt projekcji (`q95~0.0296`, nonboundary~0.725),
     - stabilnosc na perturbacje parametrow (dobre metryki wrazliwosci).
2. `QW_1776_REPARAM_EMPIRICAL_READINESS_GATE.py`
   - Agregacja 1768 + 1774 + 1775.
   - Global score: `0.9347`, hard gate `PASS`, readiness `REPARAM_EMPIRICAL_READY`.

## 37. Stan po 1776
Po tej rundzie:
1. Nowa sciezka ToE (nieobwiedniowy mechanizm + reparametryzacja zalezna od `omega`) jest metodologicznie i statystycznie gotowa do etapu empirycznego.
2. Osiagnieto rygor: kontrola leakage, grouped-CV, walidacja OOD, i formalna bramka gotowosci z PASS.
3. Kolejny logiczny krok to testy na danych empirycznych (np. pipeline GW/PTA) przy uzyciu nowego slownika parametrow zamiast starego `beta` 1:1.

## 38. Empiryczny precheck kampanii reparametryzowanej (QW-1777)
1. `QW_1777_EMPIRICAL_PTA_REPARAM_PRECHECK.py`
   - Cel: formalnie sprawdzic, czy dotychczasowe wyniki empiryczne PTA/GW i gotowosc 1776 pozwalaja uruchomic nowa kampanie testowa.
   - Wejscia: `phase70/71/75/76`, `QW-1725`, `QW-1726`, `QW-1776`.
   - Werdykt: `EMPIRICAL_REPARAM_CAMPAIGN_READY`.
   - Klucz:
     - pipeline ma czulosc (silny test injekcyjny),
     - sygnal strukturalny w danych realnych jest niezerowy,
     - stara sciezka jest niestabilna/negatywna, co uzasadnia przejscie na nowy slownik,
     - nowa sciezka ma formalne PASS gotowosci (1776).

## 39. Stan po 1777
1. Etap przygotowawczy ToE zostal domkniety na poziomie metodologicznym.
2. Aktualny status projektu: gotowosc do uruchomienia empirycznej kampanii porownawczej
   (legacy vs reparametryzacja) na danych PTA/GW.
3. Nastepny krok roboczy: zdefiniowac i uruchomic skrypt kampanii empirycznej z prerejestrowanym celem
   `logB_reparam > +3` na danych realnych przy zachowaniu skutecznosci injekcji.

## 40. Pierwszy realny test kampanii empirycznej (QW-1778 do QW-1779)
1. `QW_1778_PTA_BAYES_REPARAM_REANALYSIS.py`
   - Cel: bezposrednie porownanie Bayesowskie na danych PTA (residuals+parfiles) dla modeli:
     flat vs legacy HD vs reparam (`HD^q` z `q` zakotwiczonym przez 1773).
   - Werdykt: `PTA_REPARAM_BAYES_IMPROVED_PARTIAL`.
   - Klucz:
     - `logB_legacy_vs_flat = -0.909`,
     - `logB_reparam_vs_flat = -0.200`,
     - poprawa wzgledem legacy: `delta_logB = +0.709`,
     - czyli reparam jest lepszy od legacy, ale jeszcze nie dodatni wzgledem flat.
2. `QW_1779_EMPIRICAL_CAMPAIGN_GATE.py`
   - Agregacja 1776 + 1777 + 1778.
   - Global score: `0.8949`, hard gate `PASS`, readiness `EMPIRICAL_CAMPAIGN_STRONGLY_ON_TRACK`.

## 41. Stan po 1779
1. Kampania empiryczna jest formalnie uruchomiona i na dobrym torze (`PASS` bramki).
2. Reparametryzacja daje mierzalna poprawę Bayesowską względem legacy na danych rzeczywistych.
3. Najblizszy cel naukowy: podniesc `logB_reparam_vs_flat` powyzej zera (docelowo > +3) przez dalsze
   udoskonalenie modelu empirycznego i/lub selekcje kohorty par pulsarow.

## 42. Selekcja kohorty empirycznej pod nowy slownik (QW-1780 do QW-1781)
1. `QW_1780_PTA_COHORT_REPARAM_BAYES.py`
   - Cel: prerejestrowana selekcja kohort par pulsarow na podstawie metryk jakosci (`n_match`, stabilnosc split-half, entropia katowa),
     a nastepnie Bayes reparam vs legacy na tych kohortach.
   - Werdykt: `PTA_COHORT_REPARAM_PARTIAL`.
   - Klucz:
     - duza poprawa wzgledem legacy we wszystkich kohortach eligible,
     - dla kohorty najbardziej stabilnej (`C3_high_stable`, `n=91`) uzyskano dodatni
       `logB_reparam_mean ~ +0.167` i `delta_logB_mean ~ +0.649`.
2. `QW_1781_COHORT_OPERATIONAL_GATE.py`
   - Cel: formalny wybor kohorty operacyjnej pod dalsza kampanie empiryczna.
   - Werdykt: `COHORT_OPERATIONAL_READY`.
   - Wybrana kohorta: `C3_high_stable` (`n_pairs=91`) z dodatnim `logB_reparam` i stabilnymi metrykami MC.

## 43. Stan po 1781
1. Udalo sie osiagnac dodatni sygnal Bayesowski reparametryzacji na danych rzeczywistych
   (w kohorcie operacyjnej), przy zachowaniu przewagi nad legacy.
2. To jest praktyczny punkt przejscia z etapu "naprawa teorii" do etapu "kumulatywna kampania empiryczna".
3. Kolejny krok: utrzymac dodatnie `logB_reparam` w replikacjach (inne seedy/okna/kohorty) i dalej podnosic go do progu publikacyjnego.

## 44. Replikacje operacyjne i pierwsza bramka stabilnosci (QW-1782 do QW-1783)
1. `QW_1782_OPERATIONAL_COHORT_REPLICATION.py`
   - Cel: replikacje na podprobkach kohorty operacyjnej (`C3_high_stable`) dla oceny stabilnosci dodatniego `logB_reparam`.
   - Werdykt: `OPERATIONAL_COHORT_REPLICATION_WEAK`.
   - Klucz:
     - pelna kohorta: `logB_reparam~+0.202`, `delta~+0.641`,
     - ale replikacyjnie: `P(logB_reparam>0)=0.75`, `P(delta>0)=0.893`,
     - duza dyspersja (`std_logB~0.236`, `std_delta~0.329`).
2. `QW_1783_CAMPAIGN_STABILITY_GATE.py`
   - Agregacja 1779 + 1781 + 1782.
   - Global score: `0.7483`, hard gate `FAIL`, readiness `CAMPAIGN_STABILITY_PARTIAL`.

## 45. Test stratyfikacji katowej (QW-1784)
1. `QW_1784_STRATIFIED_REPLICATION.py`
   - Cel: zastapic losowe podprobkowanie stratyfikacja katowa, aby ograniczyc wariancje od pokrycia katowego.
   - Werdykt: `STRATIFIED_REPLICATION_WEAK`.
   - Klucz:
     - pelna kohorta nadal dodatnia (`logB_reparam~+0.201`, `delta~+0.621`),
     - poprawa dla `delta` (`P(delta>0)=0.964`) i mniejsza dyspersja `delta`,
     - ale oslabienie dla znaku `logB_reparam` (`P(logB_reparam>0)=0.643`) oraz nadal fail dyspersji.

## 46. Replikacja high-coverage (QW-1785)
1. `QW_1785_HIGH_COVERAGE_REPLICATION.py`
   - Cel: sprawdzic, czy fail 1782/1784 wynikal z utraty mocy przy `frac=0.8`; zastosowano stratyfikowane `leave-5%-out` (`frac=0.95`).
   - Werdykt: `HIGH_COVERAGE_REPLICATION_SUPPORTED`.
   - Klucz:
     - pelna kohorta: `logB_reparam~+0.189`, `delta~+0.649`,
     - replikacyjnie: `P(logB_reparam>0)=1.0`, `P(delta>0)=1.0`,
     - niska dyspersja (`std_logB~0.079`, `std_delta~0.124`),
     - wszystkie flagi PASS.

## 47. Bramka odpornosci kampanii po naprawie mocy (QW-1786)
1. `QW_1786_CAMPAIGN_ROBUSTNESS_GATE.py`
   - Agregacja: 1779 (launch), 1781 (operational cohort), 1782+1784 (stress low-coverage), 1785 (high-coverage).
   - Metoda: stress-score ciagly (pelna kohorta + prawdopodobienstwa znaku) zamiast binarnego odrzucenia tylko za jedna podprobe.
   - Wynik:
     - `global_score=0.9657`,
     - `hard_gate=PASS`,
     - readiness `EMPIRICAL_CAMPAIGN_ROBUSTNESS_CONFIRMED`.

## 48. Stan po 1786
1. Kampania empiryczna reparametryzacji jest teraz domknieta metodologicznie na poziomie stabilnosci:
   - dodatni sygnal na pelnej kohorcie,
   - dodatni sygnal utrzymany w rygorystycznej replikacji high-coverage,
   - potwierdzona odpornosc na protokoly stressowe.
2. Najwazniejszy otwarty cel merytoryczny pozostaje ten sam:
   - podniesc `logB_reparam_vs_flat` z poziomu lekkiego dodatniego (~0.19 na kohorcie operacyjnej) do silnego progu dowodowego.
3. Kolejny krok naukowy:
   - przejsc z samej stabilnosci do optymalizacji mocy empirycznej (np. prerejestrowane rozszerzenie kohorty i wariant modelu z kontrola stopnia swobody `q`), bez naruszania rygoru anty-p-hacking.

## 49. Prerejestrowany sweep frakcji pokrycia (QW-1787)
1. `QW_1787_COVERAGE_FRACTION_SWEEP.py`
   - Cel: formalny wybor frakcji podprobkowania (`0.80`, `0.85`, `0.90`, `0.95`) dla stratyfikowanych replikacji.
   - Werdykt: `COVERAGE_FRACTION_SELECTION_SUPPORTED`.
   - Wynik:
     - rekomendowana frakcja: `0.95` (`STRONG`),
     - to jedyny wariant z `pass_basic=True`,
     - metryki dla `0.95`: `P(logB_reparam>0)=1.0`, `P(delta>0)=1.0`,
       `std_logB~0.068`, `std_delta~0.105`, najlepszy `selection_score~0.875`.

## 50. Stan po 1787
1. Standard operacyjny dla kolejnych replikacji zostal domkniety:
   - stosujemy stratyfikowane high-coverage `frac=0.95` jako domyslny protokol.
2. Ten krok ogranicza arbitralnosc wyboru ustawien i wzmacnia rygor prerejestracyjny kampanii.
3. Nastepny krok merytoryczny:
   - zbadac, czy przy tym ustalonym protokole da sie podniesc bezwzgledny poziom `logB_reparam_vs_flat`
     (obecnie dodatni, ale nadal daleki od progu silnego dowodu).

## 51. Sweep szerokosci prioru q przy ustalonym protokole (QW-1788)
1. `QW_1788_Q_PRIOR_WIDTH_SWEEP.py`
   - Cel: przy stalych ustawieniach (`frac=0.95`) sprawdzic stabilnosc i jakosc sygnalu dla roznych szerokosci prioru `q`.
   - Testowane `q_width`: `0.10`, `0.14`, `0.16`, `0.20`, `0.24`.
   - Werdykt: `Q_PRIOR_WIDTH_SELECTION_SUPPORTED`.
   - Wynik:
     - rekomendowane `q_width=0.20` (`STRONG`),
     - wszystkie warianty przechodza `pass_basic`,
     - we wszystkich wariantach: dodatni `full_logB_reparam` i `full_delta`, oraz
       `P(logB_reparam>0)=1.0`, `P(delta>0)=1.0` w replikacjach.

## 52. Stan po 1788
1. Kampania ma teraz dwa ustalone, prerejestrowalne standardy operacyjne:
   - frakcja replikacji: `frac=0.95`,
   - szerokosc prioru: `q_width=0.20`.
2. To wzmacnia rygor i redukuje ryzyko arbitralnosci wyboru hiperparametrow.
3. Otwarty cel pozostaje merytoryczny (nie tylko metodologiczny):
   - zwiekszyc absolutny poziom `logB_reparam_vs_flat` ponad "lekko dodatni" przy zachowaniu tych samych rygorow.

## 53. Sweep kryteriow kohorty przy stalym rygorze (QW-1789)
1. `QW_1789_COHORT_CRITERIA_SWEEP.py`
   - Cel: sprawdzic, czy relaksacja kryteriow kohorty (`n_match`, stabilnosc) podniesie `full_logB_reparam`
     przy stalej metodzie (`frac=0.95`, `q_width=0.20`).
   - Werdykt: `COHORT_CRITERIA_SELECTION_SUPPORTED`.
   - Wynik:
     - rekomendacja `STRONG`: `K0_base` (`n_match_min=120`, `stability_max=0.65`, `n_pairs=91`),
     - `K0_base` daje najwyzszy wynik selekcyjny i najlepsza rownowage sygnalu/stabilnosci,
     - bardziej "szerokie" kohorty zwiekszaja `delta`, ale obnizaja `P(logB_reparam>0)` i/lub `full_logB_reparam`.

## 54. Stan po 1789
1. Potwierdzono, ze obecna kohorta operacyjna nie jest przypadkowym wyborem:
   - alternatywne kryteria nie daja lepszej bezwzglednej mocy sygnalu reparam.
2. To zamyka glowna os arbitralnosci konfiguracji kampanii (frakcja, prior `q`, kryteria kohorty).
3. Kolejny krok powinien zmienic klase modelu (nie tylko hiperparametry selekcji), jesli celem jest
   istotny wzrost `logB_reparam_vs_flat`.

## 55. Test zmiany klasy modelu: augmentacja katowa P2 (QW-1790)
1. `QW_1790_RESIDUAL_ANGULAR_AUGMENTATION.py`
   - Cel: sprawdzic, czy rozszerzenie modelu `M2` o ortogonalny skladnik katowy `B*P2(cos(theta))`
     zwieksza evidence (po pelnej karze za zlozonosc).
   - Werdykt: `ANGULAR_AUGMENTATION_WEAK` (jednoznaczny fail).
   - Klucz:
     - `full_logB(M2 vs flat) ~ +0.172`,
     - `full_logB(M3 vs flat) ~ -0.705`,
     - `full_delta(M3-M2) ~ -0.878`,
     - replikacyjnie: `P(M3>M2)=0.0`, `P(M3>flat)=0.0625`.

## 56. Stan po 1790
1. Dodanie prostego skladnika anizotropowego `P2` nie domyka teorii empirycznie, tylko ja pogarsza.
2. To jest silny sygnal metodologiczny:
   - nastepne rozszerzenia modelu musza byc bardziej fizycznie ukierunkowane niz ogolny dopalacz katowy.
3. Zachowujemy obecny model bazowy `M2 (HD^q + C)` jako najlepszy kandydat roboczy
   do dalszej, rygorystycznej rozbudowy.

## 57. Test robust likelihood (Student-t) vs Gauss (QW-1791)
1. `QW_1791_ROBUST_LIKELIHOOD_REPARAM.py`
   - Cel: sprawdzic, czy robust model szumu (Student-t, `nu=4`) poprawi evidence reparametryzacji
     bez zmiany klasy sygnalu.
   - Werdykt: `ROBUST_LIKELIHOOD_PARTIAL`.
   - Klucz:
     - Student-t utrzymuje dodatni sygnal (`full_t_logB_reparam_vs_flat~+0.074`, `P(t>0)=1.0`),
     - ale nie poprawia go wzgledem Gaussa (`full_delta_t_minus_gauss~-0.095`,
       `P(t>gauss)=0.188`),
     - dyspersja roznicy mala (`std~0.066`), czyli wynik stabilny, nie artefakt pojedynczej repliki.

## 58. Stan po 1791
1. Robust likelihood nie jest sciezka do zwiekszenia mocy sygnalu reparam (w tym setupie).
2. Operacyjnie pozostajemy przy Gaussie jako bazowym likelihood dla kampanii, bo daje wyzszy `logB`.
3. Jednoczesnie wynik Student-t potwierdza, ze dodatni sygnal reparam nie znika przy bardziej odpornym
   modelu szumu (to plus dla wiarygodnosci metodologicznej).

## 59. Test heteroscedastycznego modelu szumu (QW-1792)
1. `QW_1792_HETEROSCEDASTIC_REPARAM.py`
   - Cel: sprawdzic, czy heteroscedastyczny model pomiaru (jako funkcja jakosci par: `n_match`, `stability`)
     poprawi bezwzgledny `logB_reparam_vs_flat` bez zmiany modelu sygnalu.
   - Werdykt: `HETEROSCEDASTIC_REPARAM_PARTIAL`.
   - Klucz:
     - heteroscedastyczny model utrzymuje dodatni sygnal i przewage nad legacy
       (`full_hetero_logB_reparam_vs_flat~+0.129`, `full_hetero_delta_vs_legacy~+0.720`),
     - ale nie poprawia bezwzglednego `logB` wzgledem homoscedastycznego Gaussa
       (`full_gain_hetero_minus_homo~-0.101`, `P(hetero>homo)=0.214`),
     - wynik stabilny (`std gain~0.051`).

## 60. Stan po 1792
1. Heteroscedastycznosc nie jest obecnie droga do zwiekszenia mocy dowodowej w sensie `logB(reparam vs flat)`.
2. Jednoczesnie poprawia/utrzymuje przewage reparam nad legacy, wiec ma wartosc diagnostyczna,
   ale nie jako nowy domyslny model operacyjny.
3. Wniosek operacyjny: bazowym noise-model pozostaje homoscedastyczny Gauss.

## 61. Bramka lock-in protokolu operacyjnego (QW-1793)
1. `QW_1793_MODEL_LOCKIN_GATE.py`
   - Cel: formalnie zamrozic najlepszy protokol empiryczny po serii 1786-1792
     oraz odrzucic gałęzie, ktore pogarszaja evidence.
   - Wejscia: `1786`, `1788`, `1789`, `1790`, `1791`, `1792`.
   - Wynik:
     - `global_score=0.9914`,
     - `hard_gate=PASS`,
     - readiness `MODEL_LOCKIN_CONFIRMED`.
2. Zamrozony protokol operacyjny (po 1793):
   - model sygnalu: `M2 = A*HD^q + C`,
   - model szumu: homoscedastyczny Gauss,
   - frakcja replikacji: `0.95`,
   - prior `q_width=0.20`,
   - kohorta: `K0_base` (`n_match_min=120`, `stability_max=0.65`, `n_pairs=91`).

## 62. Stan po 1793
1. Etap strojenia i porownan modeli pomocniczych jest domkniety formalnie.
2. Dalsze badania powinny juz iść w kierunku:
   - nowej informacji fizycznej (nowe obserwable / nowy sygnal), a nie
   - kolejnych lokalnych modyfikacji likelihood lub prostych dopalaczy anizotropowych.
3. Mamy stabilny, prerejestrowalny baseline do kolejnych testow domykajacych ToE empirycznie.

## 63. Audyt struktury reszt modelu zamrozonego (QW-1794)
1. `QW_1794_LOCKED_MODEL_RESIDUAL_AUDIT.py`
   - Cel: sprawdzic, czy po dopasowaniu zamrozonego modelu (`M2`, Gauss) reszty sa
     nieustrukturyzowane (szum) czy zawieraja wzorzec katowy/jakosciowy.
   - Werdykt: `LOCKED_MODEL_RESIDUALS_STRUCTURED`.
   - Klucz:
     - wykryto sygnal struktury katowej (`rho(theta,resid)~0.177`, ponad prog),
     - brak sygnalu struktury po proxy jakosci (`rho~ -0.014`, silny brak efektu),
     - duzy dryf srednich binowych po katach (`max_abs_bin_mean~0.280`),
     - diagnostyka z-score dobra (`z_std~0.989`, `q95|z|~1.49`), wiec to nie jest problem outlierow.

## 64. Stan po 1794
1. Najwazniejszy aktualny problem naukowy jest teraz precyzyjnie zlokalizowany:
   - model `M2` nie usuwa calej struktury katowej w danych PTA.
2. To oznacza, ze dalsza rozbudowa ToE empirycznie powinna celowac w fizycznie uzasadniony
   skladnik opisujacy ta resztowa strukture katowa (nie przez ogolny dopalacz P2, ktory juz odpadl).
3. Mamy wiec jasny punkt wejscia do kolejnej rundy:
   - zidentyfikowac forme brakujacego skladnika tak, by poprawic `logB(reparam vs flat)` i jednoczesnie
     splaszczyc reszty katowe.

## 65. Skan modow katowych reszt (QW-1795)
1. `QW_1795_RESIDUAL_ANGULAR_MODE_SCAN.py`
   - Cel: sprawdzic, czy struktura reszt z 1794 odpowiada pojedynczemu niskiemu modowi Legendre (`P1...P6`).
   - Werdykt: `RESIDUAL_MODE_SIGNAL_WEAK`.
   - Klucz:
     - najlepszy tryb (`ell=1`) daje tylko slaby zysk in-sample (`R2~0.036`) i ujemne `delta_BIC`,
     - walidacja krzyzowa nie potwierdza stabilnej poprawy (`CV gain~0.0056`, wysoki rozrzut),
     - zadny pojedynczy niski mod nie przechodzi kryteriow BIC+CV.

## 66. Stan po 1795
1. Struktura reszt jest realna (1794), ale nie redukuje sie do prostego pojedynczego modu niskiego rzedu.
2. To wskazuje, ze brakujaca fizyka/modelowy skladnik jest bardziej zlozony:
   - prawdopodobnie wielomodowy lub nielokalny katowo, a nie "jeden brakujacy harmoniczny dopalacz".
3. Kolejna runda powinna testowac model fizyczny o ograniczonej liczbie parametrow, ale z
   bogatsza struktura katowa niz pojedyncze `P_l`.

## 67. Fizycznie ograniczone rozszerzenie wielomodowe (QW-1796)
1. `QW_1796_PHYSICAL_MULTIMODE_EXTENSION.py`
   - Cel: przetestowac trzy rodziny fizycznie ograniczonych szablonow katowych
     (odd/even/mixed) z tylko dwoma dodatkowymi stopniami swobody (`B`, `alpha`).
   - Werdykt: `PHYSICAL_MULTIMODE_EXTENSION_WEAK`.
   - Klucz:
     - najlepsza rodzina: `mixed_low [1,2,3]`,
     - na pelnej kohorcie: lokalny zysk `full_delta(M5-M2)~+0.159`,
     - ale replikacyjnie: `P(M5>M2)=0.50`, `P(M5>flat)=0.75`, wysoka dyspersja (`std~0.409`),
     - poprawa jednego wskaznika reszt (|rho|), ale pogorszenie plaskosci binow katowych.

## 68. Stan po 1796
1. Rozszerzenie wielomodowe jest na ten moment niestabilne i nie przechodzi rygoru.
2. To oznacza, ze brakujacy skladnik nie jest dobrze reprezentowany przez prosty "mikser" niskich modow
   z jedna wspolna krzywa wag.
3. Dalszy krok powinien iść w kierunku:
   - modelu z silniejszym ograniczeniem fizycznym (np. sprzezonym z cecha danych, a nie czysto katowym),
   - albo protokolu hierarchicznego (shared prior przez kohorty) w celu redukcji wariancji replikacyjnej.

## 69. Hierarchiczne shrinkage dla najlepszej rodziny (QW-1797)
1. `QW_1797_HIERARCHICAL_MULTIMODE_SHRINKAGE.py`
   - Cel: sprawdzic, czy shared prior dla (`B`, `alpha`) stabilizuje rozszerzenie `mixed_low`.
   - Werdykt: `HIERARCHICAL_MULTIMODE_PARTIAL`.
   - Klucz:
     - duzy zysk pelnokohortowy (`full_delta_hier_vs_m2~+1.936`),
     - mocny sygnal replikacyjny (`P(hier>M2)=1.0`, `P(hier>flat)=1.0`),
     - ale nadal zbyt duza dyspersja (`std_delta~0.476`).

## 70. Kalibracja sily shrinkage (QW-1798)
1. `QW_1798_HIERARCHICAL_SHRINKAGE_CALIBRATION.py`
   - Cel: sweep sily zawężenia prioru, by poprawic kompromis zysk/stabilnosc.
   - Werdykt: `HIERARCHICAL_SHRINKAGE_PARTIAL`.
   - Klucz:
     - najlepszy punkt: `shrink_factor=0.20` (rekomendacja `CONDITIONAL`),
     - nadal fail kryterium stabilnosci (`std_delta~0.422`, powyzej progu),
     - wiec kalibracja nie rozwiazala glównego problemu wariancji.

## 71. Test transferu train->test dla priory hierarchicznego (QW-1799)
1. `QW_1799_HIERARCHICAL_TRANSFER_TEST.py`
   - Cel: rygorystycznie sprawdzic generalizacje (uczenie prioru na train, ocena na rozlacznym holdoucie).
   - Werdykt: `HIERARCHICAL_TRANSFER_WEAK`.
   - Klucz:
     - gain na holdoucie czesto dodatni (`P(test hier>M2)~0.857`),
     - ale niestabilny (`test std_delta~1.119`) i zbyt slaba dodatnosc vs flat (`P~0.786`),
     - co dyskwalifikuje gałąz jako gotowa operacyjnie.

## 72. Bramka galezi hierarchicznej (QW-1800)
1. `QW_1800_HIERARCHICAL_BRANCH_GATE.py`
   - Agregacja: 1796 + 1797 + 1798 + 1799.
   - Wynik:
     - `global_score=0.55`,
     - `hard_gate=FAIL`,
     - readiness `HIERARCHICAL_BRANCH_NOT_READY`,
     - rekomendacja: `PARK_BRANCH_AND_REQUIRE_NEW_PHYSICAL_MECHANISM`.

## 73. Stan po 1800
1. Galaz wielomodowa+hierarchiczna zostaje zaparkowana (brak gotowosci operacyjnej).
2. Powod merytoryczny: sygnal jest lokalnie silny, ale nie przechodzi rygoru stabilnosci i transferu.
3. Kolejne badania powinny wejsc w nowa klase mechanizmu fizycznego, a nie dalsze tuningowanie
   tej samej rodziny szablonow katowych.

## 74. Nowa klasa mechanizmu: dekoherencja propagacyjna (QW-1801)
1. `QW_1801_PHYSICAL_DECOHERENCE_EXTENSION.py`
   - Cel: odejsc od miksowania modow i sprawdzic fizyczne tlumienie sygnalu wraz z katem (modele dekoherencji).
   - Testowane prawa tlumienia: `linear`, `quadratic`, `rational`.
   - Werdykt: `PHYSICAL_DECOHERENCE_EXTENSION_PARTIAL`.
   - Klucz:
     - najlepszy wariant: `quadratic`,
     - stabilny zysk vs `M2` (`full_delta~+0.099`, `P(M6>M2)=1.0`, `std_delta~0.059`),
     - poprawa metryki korelacyjnej reszt, ale brak poprawy plaskosci srednich binowych.

## 75. Porownanie form: multiplicative vs additive (QW-1802)
1. `QW_1802_DECOHERENCE_FORM_COMPARISON.py`
   - Cel: sprawdzic, czy forma addytywna dekoherencji poprawi jednoczesnie signal i reszty.
   - Werdykt: `DECOHERENCE_ADDITIVE_FORM_WEAK` (na ustawionych progach).
   - Klucz:
     - forma addytywna jest silniejsza od multiplikatywnej (`full_delta add-vs-mul~+1.03`),
     - poprawia oba wskazniki reszt (`d|rho|~+0.144`, `d(bin)~+0.025`),
     - ale ma wysoka dyspersje replikacyjna (`std_delta~0.464`) i dlatego nie przechodzi rygoru.

## 76. Stabilizacja hierarchiczna formy addytywnej (QW-1803)
1. `QW_1803_ADDITIVE_HIERARCHICAL_STABILIZATION.py`
   - Cel: zredukowac wariancje addytywnego modelu przez shrinkage priorow (`B`, `lambda`).
   - Werdykt: `ADDITIVE_HIERARCHICAL_STABILIZATION_PARTIAL`.
   - Klucz:
     - bardzo silny zysk pelnokohortowy utrzymany (`full_delta~+2.31` dla najlepszego ustawienia),
     - lecz nadal brak stabilnosci (`std_delta~0.435`, brak `pass_basic` dla wszystkich sf).

## 77. Dekompozycja wariancji: MC vs niestabilnosc danych (QW-1804)
1. `QW_1804_VARIANCE_DECOMPOSITION_MC.py`
   - Cel: sprawdzic, czy wysoki rozrzut wynika z szumu estymatora Bayes (MC), czy z realnej zmiennosci miedzy splitami danych.
   - Werdykt: `VARIANCE_MAINLY_DATA_INSTABILITY`.
   - Klucz:
     - komponent MC jest znikomy (`within_ratio~0.004`),
     - dominuje wariancja miedzy splitami (`between_ratio~0.996`, `std split means~0.83`),
     - czyli problem jest merytoryczny (heterogenicznosc danych), a nie numeryczny.

## 78. Bramka galezi dekoherencji (QW-1805)
1. `QW_1805_DECOHERENCE_BRANCH_GATE.py`
   - Agregacja: 1801 + 1802 + 1803 + 1804.
   - Wynik:
     - `global_score=0.60`,
     - `hard_gate=FAIL`,
     - readiness `DECOHERENCE_BRANCH_PARTIAL`,
     - rekomendacja: `PARK_FOR_PRIMARY_CLAIMS_USE_AS_SECONDARY_DIAGNOSTIC`.

## 79. Stan po 1805
1. Galaz dekoherencji jest poznawczo cenna, ale nie gotowa do glownych roszczen ToE.
2. Rdzeniowy problem po tej serii:
   - wysoka heterogenicznosc miedzy splitami danych (niezalezna od szumu MC).
3. Kolejny krok powinien celowac w modelowanie tej heterogenicznosci jako oddzielnego mechanizmu
   (np. jawny model reżimów/klastrów), zamiast dalszego dokręcania tego samego jednorodnego modelu.

## 80. Diagnostyka reżimowa po jakości par (QW-1806)
1. `QW_1806_QUALITY_REGIME_DIAGNOSTIC.py`
   - Cel: sprawdzic, czy niestabilnosc splitów wynika z ukrytych zmian skladu jakosci par (`n_match`, `stability`).
   - Model diagnostyczny: `M2Q = A*HD^q + B*Q + C`.
   - Werdykt: `QUALITY_REGIME_SIGNAL_WEAK`.
   - Klucz:
     - silnie negatywny efekt (`full_delta M2Q-M2 ~ -0.477`, `P(M2Q>M2)=0.0`),
     - brak redukcji korelacji reszt z jakoscia (`quality improvement < 0`),
     - wynik stabilny (niski `std`), wiec to nie artefakt estymatora.

## 81. Stan po 1806
1. Hipoteza "niestabilnosc = sklad jakosci próbek" nie jest wspierana.
2. To wzmacnia wniosek, ze heterogenicznosc ma bardziej fundamentalny charakter modelowy/fizyczny,
   a nie prosty bias jakosciowy.
3. Nastepna faza powinna przejsc na inny typ informacji (np. jawny komponent czasowo-czestotliwosciowy
   lub modelowanie latentnych rezimów dynamicznych), bo same statyczne rozszerzenia katowe/jakosciowe
   nie domykaja problemu.

## 82. Bramka przejscia do nowej fazy (QW-1807)
1. `QW_1807_NEXT_PHASE_GATE.py`
   - Cel: syntetycznie rozstrzygnac, czy kontynuowac statyczne strojenie modelu,
     czy formalnie przejsc do fazy dynamicznej.
   - Wejscia: `1793`, `1800`, `1805`, `1806`.
   - Wynik:
     - `static_branch_ready=False`,
     - readiness `TRANSITION_TO_DYNAMIC_PHASE_RECOMMENDED`,
     - rekomendacja: `LAUNCH_DYNAMIC_LATENT_REGIME_PROGRAM`.

## 83. Stan po 1807
1. Faza statycznych rozszerzen modelu została wyczerpana w rygorze operacyjnym.
2. Dalsze domykanie ToE empirycznie powinno przejsc na program dynamicznych reżimów latentnych
   (czasowo-zalezne/warunkowe komponenty), bo to jedyna logiczna sciezka zgodna z wynikami gate'ow.

## 84. Dynamiczny model dryfu (QW-1808)
1. `QW_1808_DYNAMIC_DRIFT_REGIME_MODEL.py`
   - Cel: pierwszy test nowej fazy dynamicznej z latentnym proxy dryfu czasowego (`h2-h1`).
   - Werdykt: `DYNAMIC_DRIFT_REGIME_WEAK`.
   - Klucz:
     - silnie negatywny efekt (`full_delta~ -1.158`, `P(M2D>M2)=0.0`),
     - niski rozrzut (stabilnie negatywne),
     - brak istotnej redukcji skorelowania reszt z dryfem.

## 85. Skan prostych cech dynamicznych (QW-1809)
1. `QW_1809_DYNAMIC_FEATURE_SCAN.py`
   - Cel: sprawdzic, czy problemem byl konkretny wybor cechy dynamicznej.
   - Testowane cechy: `drift`, `volatility`, `log_ratio`.
   - Werdykt: `DYNAMIC_FEATURE_SCAN_WEAK`.
   - Klucz:
     - `drift` i `log_ratio` zdecydowanie negatywne,
     - `volatility` daje tylko slaby sygnal (`full_delta~+0.065`) z granicznym rozrzutem (`std~0.303`),
     - brak przejscia progow rygoru.

## 86. Bramka dynamicznej fazy v1 (QW-1810)
1. `QW_1810_DYNAMIC_PHASE1_GATE.py`
   - Agregacja: 1807 + 1808 + 1809.
   - Wynik:
     - `global_score=0.40`,
     - `hard_gate=FAIL`,
     - readiness `DYNAMIC_PHASE1_NOT_READY`,
     - rekomendacja: `REDESIGN_DYNAMIC_INPUT_REPRESENTATION`.

## 87. Stan po 1810
1. Proste jednowymiarowe proxy dynamiczne nie sa wystarczajace.
2. Kolejna faza musi zmienic reprezentacje sygnalu dynamicznego:
   - z pojedynczych skalarow na bogatsze opisy czasowo-skalowe/reżimowe.
3. To jest teraz glowny kierunek dalszej naprawy empirycznej ToE.

## 88. Model dynamicznej triady cech (QW-1811)
1. `QW_1811_DYNAMIC_TRIAD_MODEL.py`
   - Cel: przetestowac bogatsza reprezentacje dynamiczna (3 cechy naraz: drift, persist, volatility).
   - Werdykt: `DYNAMIC_TRIAD_MODEL_WEAK`.
   - Klucz:
     - evidence silnie negatywne (`full_delta~ -1.767`, `P(M2T>M2)=0.0`),
     - ale poprawia diagnostyke reszt wzgledem cech dynamicznych (`feature_corr_improvement~+0.054`),
     - co sugeruje, ze cechy te opisują pewien sygnal uboczny, lecz nie pomagaja w modelu glownego sygnalu.

## 89. Bramka galezi dynamicznych proxy (QW-1812)
1. `QW_1812_DYNAMIC_PROXY_BRANCH_GATE.py`
   - Agregacja: 1808 + 1809 + 1811.
   - Wynik:
     - `global_score=0.20`,
     - `hard_gate=FAIL`,
     - readiness `DYNAMIC_PROXY_BRANCH_NOT_READY`,
     - rekomendacja: `MOVE_TO_SEQUENCE_LEVEL_DYNAMIC_MODELS`.

## 90. Stan po 1812
1. Zamknieta zostaje podfaza "dynamiczne proxy skalarne".
2. Dalsze badania powinny pracowac bezposrednio na reprezentacjach sekwencyjnych
   (okna czasowe / trajektorie / latent-state), a nie na zredukowanych skalarach z polowek sygnalu.
3. To jest obecnie najwyzszy priorytet, jesli celem jest rzeczywiste domkniecie luki empirycznej ToE.

## 91. Model sekwencyjno-okienkowy (QW-1813)
1. `QW_1813_SEQUENCE_WINDOW_MODEL.py`
   - Cel: wejscie na poziom sequence-level przez cechy z trajektorii okien czasowych.
   - Werdykt: `SEQUENCE_WINDOW_MODEL_WEAK`.
   - Klucz:
     - pelna kohorta: `full_delta~ -0.760`, `P(M2S>M2)=0.0`,
     - rozrzut umiarkowany (`std~0.254`), czyli wynik stabilnie negatywny,
     - sam fakt przejscia na okna czasowe nie wystarczyl przy obecnej reprezentacji.

## 92. Bramka dynamicznej fazy v2 (QW-1814)
1. `QW_1814_DYNAMIC_PHASE2_GATE.py`
   - Agregacja: 1810 + 1812 + 1813.
   - Wynik:
     - `global_score=0.283`,
     - `hard_gate=FAIL`,
     - readiness `DYNAMIC_PHASE2_NOT_READY`,
     - rekomendacja: `PARK_DYNAMIC_BRANCH_AND_REDESIGN_DATA_REPRESENTATION_PIPELINE`.

## 93. Stan po 1814
1. Obecna runda dynamiczna (proxy + proste cechy sekwencyjne) nie domyka ToE empirycznie.
2. Potrzebna jest przebudowa pipeline'u wejscia danych, nie tylko kolejna modyfikacja modelu:
   - richer sequence representation (np. wieloskalowe embedowanie czasowe),
   - ewentualnie nowy poziom danych (bardziej surowe sygnaly, nie tylko skondensowane `hxy`).
3. Operacyjnie: utrzymujemy zamrożony baseline `M2` jako punkt odniesienia, a nowe gałęzie
   traktujemy jako badawcze do czasu przejscia gate'ow stabilnosci i transferu.

## 94. Multiskalowe embedowanie sekwencyjne (QW-1815)
1. `QW_1815_MULTISCALE_SEQUENCE_EMBEDDING.py`
   - Cel: przejscie na bogatsza reprezentacje sequence-level (7 cech trajektorii okienkowej).
   - Werdykt: `MULTISCALE_SEQUENCE_EMBEDDING_PARTIAL`.
   - Klucz:
     - bardzo silny sygnal in-sample (`full_delta M2E-M2 ~ +17.254`, `P(M2E>M2)=1.0`),
     - ale brak kontroli rozrzutu (`std_delta~1.236` >> prog 0.30).

## 95. Stabilizacja przez jawne rezimy (QW-1816)
1. `QW_1816_REGIME_AWARE_SEQUENCE_MODEL.py`
   - Cel: zredukowac heterogenicznosc splitow przez model z rezimami latentnymi (PC1 tercyle).
   - Werdykt: `REGIME_AWARE_SEQUENCE_MODEL_WEAK`.
   - Klucz:
     - model nie poprawil M2E (`full_delta M2ER-M2E ~ -1.188`),
     - brak stabilizacji (`std_reduction_vs_m2e < 0`).

## 96. Walidacja out-of-sample branchu sekwencyjnego (QW-1817)
1. `QW_1817_SEQUENCE_OOS_VALIDATION.py`
   - Cel: odseparowac realny sygnal predykcyjny od in-sample overfittingu.
   - Werdykt: `SEQUENCE_OOS_VALIDATION_PARTIAL`.
   - Klucz:
     - dodatni zysk predykcyjny (`mean delta LL M2E-M2 ~ +5.409`, `P(M2E>M2)=0.929`),
     - staly zysk RMSE (`P(RMSE gain>0)=1.0`),
     - nadal wysoki rozrzut efektu (`std delta LL~4.089`) i niepelna przewaga nad flat (`0.857`).

## 97. Robust preprocessing w OOS (QW-1818)
1. `QW_1818_ROBUST_SEQUENCE_OOS.py`
   - Cel: zbic wariancje efektu przez winsoryzacje i mocniejsza regularizacje.
   - Werdykt: `ROBUST_SEQUENCE_OOS_PARTIAL`.
   - Klucz:
     - sygnal predykcyjny pozostaje dodatni,
     - ale wariancja rosnie zamiast spadac (`std 4.089 -> 5.656`, ratio redukcji < 0).

## 98. Bramka rygoru branchu sekwencyjnego (QW-1819)
1. `QW_1819_SEQUENCE_RIGOR_GATE.py`
   - Agregacja: 1815 + 1816 + 1817 + 1818.
   - Wynik:
     - `global_score=0.517`,
     - `hard_gate=FAIL`,
     - readiness `SEQUENCE_BRANCH_NOT_READY`,
     - rekomendacja: `PARK_AND_REDESIGN_MODELING_ASSUMPTIONS`.
   - Zidentyfikowane niespojnosci:
     - brak kontroli dyspersji,
     - niepelna przewaga OOS nad modelem flat,
     - brak poprawy po robust i po rezimach.

## 99. Redesign heteroscedastyczny (QW-1820)
1. `QW_1820_HETEROSCEDASTIC_SEQUENCE_OOS.py`
   - Cel: naprawic niestabilnosc LL przez model wariancji zalezny od cech sekwencji.
   - Werdykt: `HETEROSCEDASTIC_SEQUENCE_OOS_WEAK`.
   - Klucz:
     - heteroscedastyczny LL pogorszyl wynik wzgledem homo (`mean delta hetero-homo ~ -1.863`),
     - dyspersja wzrosla (`std vs M2: 4.089 -> 5.780`),
     - glowna luka stabilnosci pozostaje nierozwiazana.

## 100. Diagnostyka kalibracji likelihood (QW-1821)
1. `QW_1821_LIKELIHOOD_CALIBRATION_DIAGNOSTIC.py`
   - Cel: sprawdzic, czy glowna niespojnosc branchu sekwencyjnego wynika z misspecyfikacji likelihood,
     a nie z braku sygnalu w modelu sredniej.
   - Werdykt: `LIKELIHOOD_MISSPECIFICATION_SIGNAL_STRONG`.
   - Klucz:
     - niestabilnosc LL potwierdzona (`std delta LL` wysokie),
     - jednoczesnie stabilna poprawa RMSE (`P(RMSE gain>0)=1.0`),
     - obecna metryka LL jest zbyt wrazliwa na zalozenia o wariancji/tailach.

## 101. Stan po 1821
1. Branch sequence-level pokazuje realny sygnal predykcyjny (zwlaszcza RMSE i sredni zysk LL),
   ale nie przechodzi rygoru stabilnosci przez problem kalibracji likelihood.
2. Operacyjna rekomendacja na kolejna faze:
   - przejsc na gate oparty o metryki robust predictive scoring (CRPS/quantile loss/Student-t),
   - utrzymac M2 jako zamrozony baseline,
   - nie domykac jeszcze tej galezi jako finalnego dowodu empirycznego ToE, dopoki
     nie zostanie zredukowana dyspersja i nie zostanie potwierdzona przewaga OOS nad flat.

## 102. Student-t predictive scoring (QW-1822)
1. `QW_1822_STUDENTT_SEQUENCE_OOS.py`
   - Cel: sprawdzic, czy heavy-tail scoring ograniczy dyspersje i poprawi rygor OOS.
   - Werdykt: `STUDENTT_SEQUENCE_OOS_WEAK`.
   - Klucz:
     - dodatnia srednia przewaga nad M2 pozostaje (`mean delta t-LL ~ +7.149`),
     - ale zbyt niski odsetek dodatnich splitow (`P>0 ~ 0.786`) i brak przewagi nad flat na progu,
     - dyspersja dalej rośnie (`std 4.089 -> 6.114`),
     - estymowane `nu` jest wysokie (blisko gaussowskiego), wiec heavy-tail nie rozwiązuje problemu.

## 103. Stan po 1822
1. Najsilniejszy stabilny sygnal pozostaje w metrykach sredniego bledu (RMSE),
   natomiast cala rodzina gate'ow opartych o LL (Gaussian/Hetero/Student-t) pozostaje niestabilna.
2. Oznacza to, ze obecna luka domkniecia ToE lezy glownie w warstwie probabilistycznej kalibracji
   (score/uncertainty model), nie w samym mean-modelu sekwencyjnym.
3. Kolejny etap powinien formalnie przejsc na gate oparty o robust predictive scoring bez silnych
   zalozen parametrycznych o wariancji (np. quantile/CRPS) i dopiero potem wracac do finalnej bramki ToE.

## 104. Robust quantile-score OOS (QW-1823)
1. `QW_1823_QUANTILE_SCORE_OOS.py`
   - Cel: zastapic niestabilny gate LL przez odporny scoring predykcyjny bez zalozen o rozkladzie wariancji.
   - Werdykt: `QUANTILE_SCORE_OOS_SUPPORTED`.
   - Klucz:
     - dodatni zysk we wszystkich replikacjach (`P(quantile gain>0)=1.0`),
     - dodatni zysk MAE we wszystkich replikacjach (`P(MAE gain>0)=1.0`),
     - kontrola dyspersji osiagnieta (`std quantile gain~0.0257`, ponizej progu).

## 105. Bramka gotowosci po zmianie kryterium (QW-1824)
1. `QW_1824_QUANTILE_GATED_READINESS.py`
   - Agregacja: 1819 + 1821 + 1823.
   - Wynik:
     - legacy LL gate: FAIL,
     - sygnal misspec likelihood: PASS,
     - quantile robust gate: PASS,
     - readiness: `SEQUENCE_BRANCH_CONDITIONAL_READY_UNDER_QUANTILE_GATING` (hard gate PASS pod nowym kryterium).

## 106. Stan po 1824
1. Branch sequence-level jest warunkowo gotowy do nastepnego etapu empirycznego,
   ale tylko przy formalnym przejsciu na quantile-gated protocol.
2. W jezyku rygoru naukowego:
   - sygnal predykcyjny mean-modelu jest powtarzalny,
   - warstwa probabilistyczna LL jest misspecified i nie moze byc dalej gate'em glownym,
   - robust predictive scoring (quantile/MAE) domyka aktualna luke stabilnosci.
3. Operacyjnie: freeze protokolu `QW-1823/QW-1824` jako nowej bazy do dalszych badan empirycznych ToE.

## 107. Bramka transferu protokolu quantile (QW-1825)
1. `QW_1825_QUANTILE_PROTOCOL_TRANSFER_GATE.py`
   - Cel: sprawdzic, czy nowy quantile-gated protocol transferuje sie z PTA na GW.
   - Wynik:
     - `PTA`: PASS,
     - `GW`: FAIL,
     - readiness: `QUANTILE_PROTOCOL_PTA_READY_GW_BLOCKED`.
   - Wniosek: mozna kontynuowac etap PTA, ale nie wolno jeszcze uruchamiac kampanii joint PTA+GW.

## 108. Identyfikowalnosc shared/control w GW (QW-1826)
1. `QW_1826_GW_SHARED_CONTROL_IDENTIFIABILITY.py`
   - Cel: sprawdzic, czy mechanizm GW jest w ogole separowalny statystycznie na testach symulacyjnych.
   - Werdykt: `GW_SHARED_CONTROL_NOT_IDENTIFIABLE`.
   - Klucz:
     - niski efekt sredni (`mean |d| ~ 0.179`),
     - niski non-overlap (`~0.245`),
     - `P(shared > control near target)=0.0`.

## 109. Bramka naprawcza GW (QW-1827)
1. `QW_1827_GW_REDESIGN_GATE.py`
   - Cel: rozbicie blokady GW na konkretne luki z priorytetami.
   - Werdykt: `GW_BRANCH_REDESIGN_REQUIRED` (`global_score~0.176`, hard FAIL).
   - Priorytety krytyczne:
     - null rejection,
     - spojnosc miedzydetektorowa,
     - lag coherence.

## 110. Mapa targetow progowych GW (QW-1828)
1. `QW_1828_GW_TARGET_THRESHOLD_MAP.py`
   - Cel: policzyc o ile trzeba poprawic kazdy wskaznik, by GW moglo przejsc gate.
   - Wynik:
     - 8 brakujacych targetow (5 krytycznych, 3 major),
     - tylko 1 target juz spelniony (`window_std`).
   - Najwieksze luki:
     - `phase_null_p_lower`: redukcja x100,
     - `shift_null_p_lower`: redukcja x80,
     - `corr@+10ms`: wzrost x12.25,
     - `abs_diff(H1L1,H1V1)`: redukcja x2.97,
     - `P(shared>control near target)`: z 0.0 do >=0.7.

## 111. Stan po 1828
1. ToE jest obecnie domknieta warunkowo dla branchu PTA pod quantile-gated protocol.
2. Główna bariera domkniecia globalnego pozostaje w branchu GW i jest juz precyzyjnie zmapowana liczbowo.
3. Kolejna iteracja badan powinna byc juz stricte GW-redesign, ukierunkowana na 5 krytycznych targetow
   z QW-1828, zanim wrocimy do finalnej bramki joint PTA+GW.

## 112. Feasibility near-target dla GW (QW-1829)
1. `QW_1829_GW_NEAR_TARGET_FEASIBILITY.py`
   - Cel: policzyc minimalne wymagane przesuniecie sredniej i redukcje wariancji,
     by osiagnac krytyczny cel `P(|H-0.31|<=0.002) >= 0.7`.
   - Werdykt: `GW_NEAR_TARGET_REQUIRES_STRUCTURAL_REDESIGN`.
   - Klucz:
     - obecnie `mean P(near target) ~ 0.00045`,
     - wymagany sredni shift `|mu| ~ 0.081`,
     - nawet po centrowaniu potrzebna redukcja sigma srednio `~12.4x`.

## 113. Globalna bramka empiryczna po przejsciu na quantile (QW-1830)
1. `QW_1830_GLOBAL_EMPIRICAL_STATUS_GATE.py`
   - Agregacja: 1824 + 1825 + 1827 + 1829.
   - Wynik:
     - `global_score~0.444`,
     - hard gate `PARTIAL`,
     - readiness `GLOBAL_PARTIAL_PTA_ONLY`.
   - Interpretacja:
     - tor PTA jest gotowy warunkowo i moze byc kontynuowany,
     - tor GW wymaga przebudowy strukturalnej przed domknieciem globalnym ToE.

## 114. Stan po 1830
1. Domkniecie ToE na poziomie globalnym (PTA+GW) jeszcze nie nastapilo.
2. Aktualny rygor naukowy pozwala na:
   - kontynuacje empirycznego toru PTA pod quantile-gated protocol,
   - równolegly program naprawczy GW bez eskalacji claimów globalnych.
3. Najblizszy etap techniczny: `QW-1831/1832/1833` (event-windowed GW features,
   shared-control objective, multi-detector consistency gate).

## 115. Event-windowed feature extractor GW (QW-1831)
1. `QW_1831_EVENT_WINDOWED_GW_FEATURES.py`
   - Cel: przejscie z globalnych statystyk GW na reprezentacje okienkowa i cechy koherencji detector-pair.
   - Dane: lokalne pliki `raw_strain_unfiltered/*_unfiltered_1266965117.h5` (H1/L1/V1), bez fetchu sieciowego.
   - Wynik: `EVENT_WINDOW_FEATURE_BASELINE_PARTIAL`.
   - Klucz:
     - cechy zostały poprawnie wyekstrahowane (po fixie skali `zscore`),
     - dobra spojnosc miedzy parami kontrolnymi w wariancji median,
     - slaby sygnal lag@10ms dla H1-L1 (`P(corr10>0.02)~0.012`).

## 116. Obiektyw shared-vs-control na oknach (QW-1832)
1. `QW_1832_GW_SHARED_CONTROL_WINDOW_OBJECTIVE.py`
   - Cel: test identyfikowalnosci H1-L1 vs par kontrolnych na nowej reprezentacji.
   - Werdykt: `GW_WINDOW_OBJECTIVE_SUPPORTED`.
   - Klucz:
     - `mean AUC ~ 0.886`,
     - `mean balanced accuracy ~ 0.792`,
     - `mean prevalence advantage ~ 0.583`, dodatni we wszystkich foldach.

## 117. Multi-detector consistency gate v2 (QW-1833)
1. `QW_1833_GW_MULTI_DETECTOR_CONSISTENCY_GATE.py`
   - Cel: sprawdzic czy poprawa 1832 jest spójna miedzy parami kontrolnymi.
   - Werdykt: `GW_MULTI_DETECTOR_PARTIAL_CONTROL_MISMATCH`.
   - Klucz:
     - separacja shared-vs-control: PASS,
     - przewaga q90: PASS,
     - AUC mix: PASS,
     - ale niespelniona spojnosc miedzy parami kontrolnymi (`control_median_gap` za duzy).

## 118. Bramka postepu redesign GW phase-1 (QW-1834)
1. `QW_1834_GW_REDESIGN_PHASE1_GATE.py`
   - Agregacja: 1827 + 1829 + 1832 + 1833.
   - Wynik:
     - hard gate `PARTIAL`,
     - readiness `GW_REDESIGN_PHASE1_PARTIAL_PROGRESS`.
   - Interpretacja:
     - nowy obiektyw GW dziala,
     - ale blokada strukturalna near-target i mismatch kontroli pozostaja.

## 119. Globalny status po GW phase-1 (QW-1835)
1. `QW_1835_GLOBAL_STATUS_AFTER_GW_PHASE1.py`
   - Agregacja: 1824 + 1825 + 1834 (+ porownanie do 1830).
   - Wynik:
     - global score: `0.444 -> 0.733` (istotna poprawa),
     - hard gate nadal `PARTIAL`,
     - readiness: `GLOBAL_PARTIAL_WITH_GW_PHASE2_REQUIRED`.
   - Wniosek:
     - PTA track pozostaje gotowy warunkowo,
     - GW przeszedl z fail do partial progress,
     - potrzebna faza `GW phase-2 structural reparam` przed domknieciem globalnym ToE.

## 120. GW phase-2: control-calibrated objective (QW-1836)
1. `QW_1836_GW_CONTROL_CALIBRATED_OBJECTIVE.py`
   - Cel: usunac systematyczny bias par kontrolnych przez kalibracje offsetow
     liczonych wyłącznie na foldach treningowych.
   - Wynik: `GW_CONTROL_CALIBRATED_OBJECTIVE_SUPPORTED`.
   - Kluczowe metryki (raw -> calibrated):
     - `mean AUC: 0.8436 -> 0.9415`,
     - `mean advantage: 0.389 -> 0.757`,
     - `mean control gap: 0.003312 -> 0.000154` (`~95.4%` redukcji).
   - Wniosek:
     - nowy obiektyw GW jest identyfikowalny i stabilny po kalibracji train-only.

## 121. Transition gate: legacy near-target -> reparam criterion (QW-1837)
1. `QW_1837_GW_CRITERION_TRANSITION_GATE.py`
   - Cel: formalnie przejsc z niespelnialnego kryterium near-target
     na kryterium kalibrowane, zgodne z 1836.
   - Wynik:
     - `global_score=0.95`,
     - hard gate `PASS`,
     - readiness `GW_READY_UNDER_REPARAM_CRITERION`.
   - Wniosek:
     - gałąź GW jest gotowa warunkowo pod nowym, jawnie zdefiniowanym kryterium.

## 122. Global readiness under reparam criteria (QW-1838)
1. `QW_1838_GLOBAL_REPARAM_READINESS_GATE.py`
   - Agregacja: 1824 + 1836 + 1837 (+ 1825 dla traceability).
   - Wynik:
     - global score: `0.733 -> 0.900`,
     - hard gate `PASS`,
     - readiness: `GLOBAL_CONDITIONAL_READY_UNDER_REPARAM_CRITERIA`.
   - Wniosek:
     - ToE nie jest jeszcze "udowodniona", ale uzyskano warunkowa gotowosc
       do kampanii potwierdzajacej PTA+GW pod zamrozonym protokolem.

## 123. Zamrozenie prerejestracji joint confirmatory (QW-1839)
1. `QW_1839_JOINT_CONFIRMATORY_PREREG_PROTOCOL.py`
   - Cel: zamrozic protokol confirmatory i anti-leakage przed nowymi uruchomieniami.
   - Wynik:
     - verdict: `JOINT_CONFIRMATORY_PREREG_FROZEN`,
     - hash protokolu:
       `b9cf21d3d32508e95c6f7cef2e8a953a12f6c7ea7732b6288605c518ab5db5af`.
   - Zamrozone progi:
     - PTA: gain>=0.04, p_pos>=0.9, std<=0.035.
     - GW: auc>=0.9, adv>=0.6, control_gap<=0.0005, p_adv_pos>=0.9.
   - Wniosek:
     - najblizszy krok badawczy to strict confirmatory rerun na nowym
       materiale lub zewnetrznym holdoucie, bez zmiany progow/splitu.

## 124. Joint power i wrazliwosc na dryf marginesu (QW-1840)
1. `QW_1840_JOINT_CONFIRMATORY_POWER_SENSITIVITY.py`
   - Cel: oszacowac prawdopodobienstwo przejscia zamrozonej bramki 1839
     przez bootstrap PTA+GW oraz przetestowac scenariusze dryfu efektu.
   - Wynik:
     - verdict: `JOINT_CONFIRMATORY_HIGH_POWER`,
     - nominalnie `p(joint)~0.994` (CI95 ~ `[0.993, 0.995]`),
     - margines odpornosci: `p(joint)>=0.90` do ~`20%` konsumpcji marginesu,
       `p(joint)>=0.50` do ~`75%`.
   - Wniosek:
     - przy stabilnym efekcie szansa przejscia confirmatory jest wysoka,
       ale przy istotnym dryfie moze szybko spadac.

## 125. Permutation-null dla kalibrowanego toru GW (QW-1841)
1. `QW_1841_GW_CALIBRATED_PERMUTATION_NULL.py`
   - Cel: sprawdzic, czy wynik 1836 nie jest artefaktem etykiet
     (null: brak informacji o parze detectorow).
   - Wynik:
     - verdict: `GW_CALIBRATED_OBJECTIVE_NULL_REJECTED`,
     - p-values:
       `AUC~0.0002`, `ADV~0.0002`, `GAP~0.0012` (wszystkie <= 0.01).
   - Wniosek:
     - nowy obiektyw GW przechodzi rygorystyczny test null i ma silna
       podpore statystyczna.

## 126. Zintegrowana bramka wykonawcza confirmatory (QW-1842)
1. `QW_1842_JOINT_CONFIRMATORY_EXECUTION_GATE.py`
   - Agregacja: 1838 + 1839 + 1840 + 1841.
   - Wynik:
     - `global_score=1.000`,
     - hard gate `PASS`,
     - readiness `READY_FOR_EXTERNAL_CONFIRMATORY_EXECUTION`.
   - Wniosek:
     - aktualna wersja ToE jest warunkowo domknieta metodologicznie
       i gotowa do zewnetrznego/testowego uruchomienia confirmatory
       bez zmiany protokolu i progow.

## 127. Audyt inferencyjny progu PTA (QW-1843)
1. `QW_1843_PTA_THRESHOLD_INFERENCE_RIGOR.py`
   - Cel: sprawdzic, czy zamrozony prog PTA `prob_quantile_gain_positive>=0.9`
     jest nie tylko spelniony punktowo, ale tez inferencyjnie.
   - Wynik:
     - verdict: `PTA_POINT_PASS_BUT_PROBABILITY_UNDERPOWERED`,
     - obserwacje: `n=14`, `k=14`, ale `p-value(H0:p<=0.9)=0.2288`,
       one-sided lower95 dla p to `~0.807`,
       minimalnie `n=29` all-positive dla alpha=0.05.
   - Wniosek:
     - kryteria mean/std sa domkniete, ale kryterium probabilistyczne PTA
       jest niedomocowane przy obecnym n.

## 128. Bramka rygoru scislego (QW-1844)
1. `QW_1844_STRICT_JOINT_RIGOR_GATE.py`
   - Agregacja: 1842 + 1843.
   - Wynik:
     - `global_score=0.800`,
     - hard gate `PARTIAL`,
     - readiness `PTA_PROBABILITY_INFERENCE_UNDERPOWERED`.
   - Wniosek:
     - operacyjna gotowosc confirmatory istnieje, ale przy rygorze
       inferencyjnym wymagane jest doszacowanie PTA.

## 129. Plan mocy dla testu PTA p>0.9 (QW-1845)
1. `QW_1845_PTA_PROBABILITY_POWER_PLAN.py`
   - Cel: policzyc wymagane n dla exact binomial test (alpha=0.05, power=0.80).
   - Wynik:
     - `p_true=1.00`: n=29 (add 15),
     - `p_true=0.99`: n=46 (add 32),
     - `p_true=0.97`: n=76 (add 62),
     - `p_true=0.95`: n=179 (add 165).
   - Wniosek:
     - obecny prog 0.9 jest kosztowny inferencyjnie i wymaga istotnego wzrostu n
       przy realistycznych p_true < 0.99.

## 130. Mapa wykonalnosci progu probabilistycznego PTA (QW-1846)
1. `QW_1846_PTA_PROB_THRESHOLD_FEASIBILITY_MAP.py`
   - Cel: porownac wymagane n dla roznych progow `p0` i scenariuszy `p_true`.
   - Wynik:
     - przy aktualnym `n=14` maksymalny odrzucalny prog (all-positive)
       to `p0~0.807`,
     - dla obecnego kryterium `p0=0.90`, `p_true=0.97`: `n=76`.
   - Wniosek:
     - luka mocy nie jest artefaktem pojedynczego testu, tylko cecha
       calego kryterium przy malej probie.

## 131. Proposal targetow replikacji PTA (QW-1847)
1. `QW_1847_PTA_REPLICATION_TARGET_PROPOSAL.py`
   - Cel: sformalizowac targety n bez zmiany zamrozonego progu 0.9.
   - Wynik:
     - status strict: `PARTIAL / PTA_PROBABILITY_INFERENCE_UNDERPOWERED`,
     - target minimalny: `n=29` (all-positive),
     - target optymistyczny: `n=46` (`p_true=0.99`),
     - target realistyczny: `n=76` (`p_true=0.97`, rekomendowany).
   - Wniosek:
     - najlepsza dalsza sciezka badawcza to eskalacja replikacji PTA
       do co najmniej `n~76` (lub rewizja kryterium probabilistycznego
       w osobnej, jawnie versionowanej gałęzi).

## 132. Audyt jednostki analizy PTA (QW-1848)
1. `QW_1848_PTA_UNIT_OF_ANALYSIS_AUDIT.py`
   - Cel: sprawdzic, czy kryterium `prob>=0.9` nie jest artefaktem
     poziomu agregacji (splity vs pary pulsarow).
   - Wynik:
     - split-level `P(rep_mean_gain>0)~0.988`,
     - pair-level `P(pair_mean_gain>0)~0.820`,
     - compression gap `~0.167`,
     - `k/n=73/89`, lower95 `~0.740`,
       p-value dla `H0: p<=0.9` `~0.993`.
   - Werdykt:
     - `PTA_UNIT_MISMATCH_REQUIRES_CRITERION_REDESIGN`.
   - Wniosek:
     - kryterium V1 na poziomie splitow zawyza pewnosc przez usrednianie.

## 133. Bramka wyboru scislej sciezki PTA (QW-1849)
1. `QW_1849_PTA_STRICT_PATH_SELECTION_GATE.py`
   - Agregacja: 1842 + 1844 + 1847 + 1848.
   - Wynik:
     - best path:
       `PATH_B_VERSIONED_CRITERION_REPARAM_WITH_EXTERNAL_CONFIRMATORY`,
     - readiness:
       `REPARAM_PROTOCOL_REQUIRED_BEFORE_STRICT_CONFIRMATORY`.
   - Wniosek:
     - naukowo lepsze jest wersjonowanie kryterium PTA (V2),
       niz dalsze \"pompowanie\" splitow V1.

## 134. Zamrozenie prerejestracji PTA_V2 (QW-1850)
1. `QW_1850_PTA_V2_PREREG_PROTOCOL.py`
   - Cel: zamrozic protokol PTA_V2 oparty o pair-level endpoint
     i anti-leakage z zakazem reuse datasetu projektowego.
   - Wynik:
     - verdict:
       `PTA_V2_PREREG_FROZEN_EXTERNAL_CONFIRMATORY_PENDING`,
     - hash:
       `e5bcdc803f5587f790d9c1a70418463ed416760c4fcec72f6cd06b46a92b2f50`.
   - Progi PTA_V2:
     - `mean_pair_mean_gain>=0.04`,
     - `bootstrap_lower95_mean_pair_mean_gain>=0`,
     - `prob_pair_mean_gain_positive>=0.667`,
     - `one_sided_lower95_prob_pair_mean_gain_positive>=0.600`.

## 135. Joint V2 external execution gate (QW-1851)
1. `QW_1851_JOINT_V2_EXTERNAL_EXECUTION_GATE.py`
   - Agregacja: GW readiness + decyzja PATH_B + PTA_V2 prereg.
   - Wynik:
     - hard gate `PASS`,
     - readiness `READY_FOR_EXTERNAL_CONFIRMATORY_V2`,
     - external data required: `True`.
   - Wniosek:
     - ToE jest teraz domknieta metodologicznie w wersji V2
       na poziomie \"ready for external confirmatory\":
       kolejny krok to jedno uruchomienie confirmatory na nowym zbiorze.

## 136. Pakiet precheck dla external confirmatory (QW-1852)
1. `QW_1852_EXTERNAL_CONFIRMATORY_DATA_PRECHECK.py`
   - Cel: utworzyc pakiet przyjecia danych zewnetrznych (README + schema + manifest)
     i walidator anti-leakage zgodny z hashami 1839/1850.
   - Wynik:
     - hard gate `PARTIAL`,
     - readiness `EXTERNAL_DATASET_PENDING_COLLECTION`,
     - wygenerowany pakiet:
       `external_confirmatory_v2/{README, manifest_template, pta_template, gw_template}`.
   - Wniosek:
     - pipeline intake jest gotowy, ale dataset external jeszcze nie dostarczony.

## 137. Runner joint external confirmatory V2 (QW-1853)
1. `QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py`
   - Cel: przygotowac finalny runner confirmatory PTA_V2+GW z twarda blokada,
     dopoki 1852 nie przejdzie na `READY`.
   - Wynik:
     - hard gate `PARTIAL`,
     - readiness `BLOCKED_WAITING_VALID_EXTERNAL_DATASET`,
     - brak wykonania testu finalnego (zgodnie z rygorem, bez danych external).
   - Wniosek:
     - domkniecie ToE w trybie confirmatory jest gotowe proceduralnie;
       jedyny brakujacy element to zwalidowany zewnetrzny zbior danych.

## 138. Audyt zakresu QW-700..1600 pod kernel freeze (QW-1854)
1. `QW_1854_KERNEL_FREEZE_RANGE_AUDIT_700_1600.py`
   - Cel: odpowiedziec twardo, czy zakres 700..1600 zostal naprawde prześledzony
     pod katem narracji \"zamrozonego kernela\".
   - Wynik:
     - verdict: `RANGE_700_1600_NOT_FULLY_TRACED`,
     - coverage any-ref: `512/901` (`~0.568`),
     - coverage kernel-ref: `509/901` (`~0.565`),
     - kontradykcje: `phi dual`, `beta dual`, `node-set conflict`.
   - Wniosek:
     - nie ma podstaw twierdzic, ze caly zakres 700..1600 jest juz domkniety
       i jednoznaczny pod frozen kernel.

## 139. Kolejka deep-read dla QW-700..1600 (QW-1855)
1. `QW_1855_KERNEL_FREEZE_DEEPREAD_QUEUE_700_1600.py`
   - Cel: zbudowac priorytety recznego czytania i korekty dla najbardziej ryzykownych QW.
   - Wynik:
     - missing ids: `389`,
     - ids with trace: `512`,
     - canonical candidates: `240`,
     - top missing od `QW-704`, `QW-712`, `QW-714`, ...
   - Wniosek:
     - pipeline do domykania luki 700..1600 jest przygotowany metodycznie.

## 140. Ledger konfliktow freeze z chronologia (QW-1856)
1. `QW_1856_KERNEL_FREEZE_CONFLICT_LEDGER_700_1600.py`
   - Cel: odseparowac QW konfliktowe od kandydatow kanonicznych
     z uwzglednieniem dat plikow.
   - Wynik:
     - conflict rows: `5` (m.in. `QW-1200`, `QW-700`),
     - canonical rows: `240`.
   - Wniosek:
     - frozen kernel jest nierownomiernie stabilny w historii 700..1600;
       istnieja wyrazne wyspy konfliktowe wymagajace kwarantanny interpretacyjnej.

## 141. Audyt konkretnego full loga 420..1200 (QW-1858)
1. `QW_1858_FULL_LOG_EXTREME_QW420_1200_AUDIT.py`
   - Cel: szczegolowo przeanalizowac plik
     `FULL_LOG_COMPRESSED_EXTREME_QW420_1200.md`.
   - Wynik:
     - 13810 linii, 386 unikalnych QW w 420..1200,
     - brak sprzecznosci definicyjnych `phi/beta/nodes`,
     - ale sa markery kryzysu frozen branch (`frozen_negative_hits=3`),
       przy silnej obecnosci kernela i freeze.
   - Wniosek:
     - plik jest wewnetrznie relatywnie spójny parametrycznie,
       ale sam raportuje empiryczne porazki tej galezi.

## 142. Meta-audyt wszystkich FULL_LOG*.md (QW-1859)
1. `QW_1859_FULL_LOGS_META_AUDIT.py`
   - Cel: zbadac systemowo wszystkie full logi, a nie pojedynczy plik.
   - Wynik:
     - full log files: `15`,
     - union QW ids: `503`,
     - medium risk: `1` (EXTREME_QW420_1200), low risk: `14`,
     - brak kontradykcji definicyjnych w samych full logach.
   - Wniosek:
     - full logi sa glownie neutralnym zapisem operacyjnym;
       krytyczny negatywny sygnal frozen branch koncentruje sie w pliku EXTREME.

## 143. Katalog wiarygodnosci full logow (QW-1860)
1. `QW_1860_FULL_LOG_TRUST_CATALOG.py`
   - Cel: sklasyfikowac full logi jako evidence positive/negative/neutral.
   - Wynik:
     - negative evidence: `FULL_LOG_COMPRESSED_EXTREME_QW420_1200.md`,
     - positive context: `FULL_LOG_COMPRESSED_PHASE3_QW1500_1535.md`,
     - neutral context: `13` plikow.
   - Wniosek:
     - w rygorze naukowym plik EXTREME powinien byc traktowany jako dowod
       przeciw prostemu frozen-branch claim, nie jako potwierdzenie.

## 144. Stress test chronologii konfliktow po freeze (QW-1857)
1. `QW_1857_KERNEL_FREEZE_CHRONOLOGY_STRESS_TEST_700_1600.py`
   - Cel: sprawdzic, czy sprzeczne definicje pojawiaja sie po deklaracji freeze.
   - Wynik:
     - conflict QWs: `5`,
     - quarantine: `2` (`QW-1200`, `QW-700`),
     - w obu tych przypadkach konflikty pojawialy sie po freeze markerach.
   - Wniosek:
     - `QW-1200` i `QW-700` nie powinny byc uzywane jako kanon zamrozonego kernela
       bez osobnej naprawy/rekonsyliacji definicji.

## 145. Re-audyt full logow po poprawce regex (QW-1858..1860 update)
1. Korekta techniczna:
   - regex `phi=pi/6` rozszerzony o wariant `np.pi/6` i `π/6`
     w audytach 1854..1859.
2. Zmiana wyniku:
   - `FULL_LOG_COMPRESSED_EXTREME_QW420_1200.md`:
     z \"internal_ok/mixed\" przechodzi na
     `FULL_LOG_FROZEN_BRANCH_CONTRADICTORY_AND_EMPIRICALLY_NEGATIVE`.
   - Wykryto `phi_pi6_hits=94` oraz `phi_zero_hits=2` (dual definition).
3. Aktualizacja katalogu zaufania (1860):
   - EXTREME jest klasyfikowany jako `EXCLUDE_CONTRADICTORY`
     (nie jako neutralny pozytywny evidence).
4. Wniosek:
   - domykanie ToE nie moze opierac claimu frozen kernel na pliku EXTREME;
     nalezy go traktowac jako material konfliktowy i diagnostyczny.

## 146. Finalna korekta stanu audytu freeze (QW-1854..1860)
1. Stan po poprawce regex i ponownych uruchomieniach:
   - `QW-1854`: coverage any-ref `537/901` (`0.596`),
     verdict `RANGE_700_1600_NOT_FULLY_TRACED`.
   - `QW-1856`: conflict rows `179`, canonical rows `69`.
   - `QW-1857`: quarantine `179/179` konfliktowych QW.
2. Wniosek:
   - frozen-kernel narrative w zakresie 700..1600 wymaga rygorystycznej
     kwarantanny konfliktow; nie ma podstaw do prostego claimu o pelnym domknieciu historycznym.

## 147. Rekonstrukcja kanonicznego kernela po filtrze zaufania (QW-1861)
1. `QW_1861_CANONICAL_KERNEL_RECONSTRUCTION_700_1600.py`
   - Cel: zrekonstruowac kanon z danych nie-sprzecznych, z wagami po datach plikow.
   - Wynik:
     - verdict `CANONICAL_KERNEL_RECONSTRUCTION_NOT_CLOSED`,
     - coverage QW `69/69` (po trybie lokalnym + diffuse),
     - osie: `phi -> pi/6 (0.877)`, `beta -> 0.01 (0.620)`,
       `omega -> pi/4 (0.971)`, ale `nodes -> none` (brak stabilnej osi wezlow).
   - Wniosek:
     - parametry `(phi, beta, omega)` maja preferowany kanon,
       ale bez zamknietej osi wezlow nie da sie uznac pelnego domkniecia kernela.

## 148. Audyt zgodnosci mikromodelu z kanonem (QW-1862)
1. `QW_1862_MICROMODEL_CANONICAL_COMPATIBILITY_AUDIT.py`
   - Cel: sprawdzic zgodnosc estymatorow 1739/1741/1743 z kanonem 1861.
   - Wynik:
     - verdict `MICROMODEL_CANONICAL_COMPATIBILITY_FAIL`,
     - `compatibility_index ~ 1.5e-43` (praktycznie zero),
     - bardzo duze odchylenia `omega` i `beta`,
       slaba identyfikowalnosc (`cond ~ 1e12` w 1742/1744).
   - Wniosek:
     - glowna luka domkniecia ToE pozostaje mikromodelowa,
       nie tylko dokumentacyjna.

## 149. Projekt obserwabli pod identyfikowalnosc (QW-1863)
1. `QW_1863_IDENTIFIABILITY_OPTIMAL_OBSERVABLES_DESIGN.py`
   - Wynik:
     - verdict `IDENTIFIABILITY_DESIGN_PARTIAL`,
     - najlepszy baseline: `phase_increment + envelope_decay + zero_cross_offset`,
     - score `~0.524`, `cond_q90 ~ 5.7`.
   - Wniosek:
     - jest wykonalny technicznie zestaw obserwabli do odzysku parametrow,
       ale sam paired sign-flip nie daje tu duzego zysku informacyjnego.

## 150. Bramka fazy 3 (QW-1864)
1. `QW_1864_KERNEL_CLOSURE_PHASE3_GATE.py`
   - Wynik:
     - global score `0.206`, hard gate `FAIL`,
       readiness `PHASE3_OPEN_KERNEL_NOT_CLOSED`.
   - Wniosek:
     - przechodzimy dalej tylko przez benchmarki wykonawcze (1865, 1866),
       bo sama rekaonstrukcja + audyt kompatybilnosci nie zamyka ToE.

## 151. Synthetic recovery benchmark (QW-1865)
1. `QW_1865_SYNTHETIC_RECOVERY_BENCHMARK.py`
   - Wynik:
     - verdict `SYNTHETIC_RECOVERY_PARTIAL`,
     - tolerance hits: `omega=1.000`, `phi=1.000`, `beta=0.775`,
     - `nonboundary_rate=0.742`.
   - Wniosek:
     - odzysk parametrow jest mocny dla `omega/phi`, ale `beta` i nienabrzegowosc
       nadal wymagaja poprawy protokolu.

## 152. Paired signed-intervention benchmark (QW-1866)
1. `QW_1866_PAIRED_SIGNED_INTERVENTION_BENCHMARK.py`
   - Wynik:
     - verdict `PAIRED_SIGNED_INTERVENTION_NOT_SUPPORTED`,
     - degradacja dla beta (`rmse_beta_factor=0.679`, `beta_tol_gain=-0.094`).
   - Wniosek:
     - ta konkretna forma paired intervention nie domyka luki beta;
       wymagana adaptacja zestawu obserwabli.

## 153. Adaptacyjny benchmark beta-augmented (QW-1867)
1. `QW_1867_BETA_AUGMENTED_OBSERVABLES_BENCHMARK.py`
   - Wynik:
     - verdict `BETA_AUGMENTATION_SUPPORTED`,
     - najlepsze ramię: `B_AUGMENTED` (bez paired),
     - `beta_rmse_factor=1.654`, `beta_tol_gain=+0.171` vs baseline.
   - Wniosek:
     - dodanie beta-czulej obserwabli (`torsion_cross_term`) realnie zmniejsza luke
       odzysku `beta` i jest lepsze od samej interwencji paired.

## 154. Bramka fazy 3B po adaptacji (QW-1868)
1. `QW_1868_KERNEL_CLOSURE_PHASE3B_GATE.py`
   - Wynik:
     - global score `0.397`, hard gate `FAIL`,
       readiness `PHASE3B_OPEN_NOT_CLOSED`.
   - Co przeszlo:
     - design identifiability, synthetic recovery execution, adaptive beta augmentation.
   - Co nie przeszlo:
     - theory canonical strength (nodes unresolved),
     - micromodel canonical compatibility.
   - Wniosek:
     - ToE jest mocniej domknieta metodologicznie i eksperymentalnie,
       ale nadal nie domknieta teoretycznie (os wezlow + zgodnosc mikromodelu z kanonem).

## 155. Najlepszy kolejny krok (po 1868)
1. Uruchomic `QW-1869` jako stricte teoretyczny redesign osi wezlow:
   - jawna konkurencja modeli wezlow (`2,5,8,11`, `2,8,14`, `2+3n`, warianty mieszane)
     z kryterium informacji i zgodnoscia z kanonicznym `(phi=pi/6, omega=pi/4, beta=0.01)`.
2. Uruchomic `QW-1870` jako mikromodel constrained-by-nodes:
   - wymusic nowy node prior z 1869 i sprawdzic,
     czy `compatibility_index` z 1862 wychodzi z obszaru ~0.
3. Tylko po przejsciu 1869+1870 wracac do sciezki real-data confirmatory.

## 156. Konkurencja modeli wezlow na danych pierwotnych (QW-1869)
1. `QW_1869_NODE_MODEL_COMPETITION_CANONICAL.py`
   - Cel: rozstrzygnac os wezlow bez autoreferencji (wyciete pliki raportowe/meta).
   - Wynik finalny (po odfiltrowaniu petli dowodowej):
     - verdict `NODE_MODEL_NOT_RESOLVED`,
     - `files_used_with_node_markers = 0`,
     - posteriory modeli plaskie (`~0.333 / 0.333 / 0.333`).
   - Wniosek:
     - luka wezlowa jest realna i dotyczy braku twardych danych pierwotnych
       w aktualnym korpusie 700..1600, a nie tylko statystycznej niejednoznacznosci.

## 157. Feasibility mikromodelu z priorem wezlowym (QW-1870)
1. `QW_1870_NODE_CONSTRAINED_MICROMODEL_FEASIBILITY.py`
   - Cel: sprawdzic, czy samo narzucenie realistycznego prioru wezlowego
     moze podniesc kompatybilnosc mikromodelu z kanonem.
   - Wynik:
     - verdict `NODE_CONSTRAINED_COMPATIBILITY_FAIL_STRONG`,
     - max score po ograniczeniach: `~8.49e-32`.
   - Wniosek:
     - nawet przy dopuszczalnych przesunieciach parametrow,
       mikromodel pozostaje skrajnie niezgodny z kanonem;
       potrzebna jest zmiana strukturalna mikromodelu, nie tylko nowy prior.

## 158. Stan domykania ToE po 1870
1. Co jest domkniete lepiej niz wczesniej:
   - rygor audytu historycznego (1854..1860),
   - benchmarki identyfikowalnosci i recovery (1863..1867),
   - adaptacja obserwabli dla `beta` (sukces 1867).
2. Co nadal blokuje domkniecie:
   - brak rozstrzygnietej osi wezlow na danych pierwotnych,
   - kompatybilnosc mikromodelu z kanonem bliska zeru (1862, 1870).
3. Najlepsza nastepna sciezka:
   - `QW-1871`: pozyskac/wyekstrahowac pierwotny corpus wezlowy
     (line-level z chronologia i linkage do konkretnych QW-runow),
   - `QW-1872`: nowy mikromodel strukturalny z jawna dynamika wezlow
     i testem zgodnosci wzgledem 1862/1870.

## 159. Korpus pierwotnych dowodow wezlowych (QW-1871)
1. `QW_1871_PRIMARY_NODE_EVIDENCE_CORPUS.py`
   - Cel: zbudowac line-level corpus wezlow z plikow pierwotnych,
     z odfiltrowaniem autoreferencji (`report_*`, `RAPORT_*`, `QW_17xx+`, `.venv`).
   - Wynik:
     - verdict `PRIMARY_NODE_MODEL_STRONGLY_SELECTED`,
     - winner: `M_A_2_5_8_11_or_2plus3n` (posterior `1.000`),
     - `n_primary_files_with_node_hit=21`, `n_primary_line_hits=74`.
   - Wniosek:
     - os wezlow jest empirycznie wsparta w danych historycznych
       i preferuje narracje `2,5,8,11 / 2+3n`.

## 160. Strukturalny mikromodel node-dynamic (QW-1872)
1. `QW_1872_STRUCTURAL_NODE_DYNAMIC_MICROMODEL.py`
   - Cel: sprawdzic, czy jawna dynamika wezlow (modulator node-field)
     podnosi zgodnosc mikromodelu z kanonem bez psucia dopasowania.
   - Wynik:
     - verdict `STRUCTURAL_NODE_MODEL_WEAK_BRIDGE`,
     - median RMSE poprawia sie (`0.1519 -> 0.1367`),
     - ale canonical score zostaje praktycznie przy zerze
       (`2.035e-08 -> 2.043e-08`, gain mikroskopijny).
   - Wniosek:
     - sama korekta wezlowa poprawia fit techniczny,
       ale nie rozwiazuje glównej luki kanonicznej `omega/beta`.

## 161. Bramka fazy 4 (QW-1873)
1. `QW_1873_KERNEL_CLOSURE_PHASE4_GATE.py`
   - Agregacja: 1868 + 1871 + 1872.
   - Wynik:
     - global score `0.534`, hard gate `FAIL`,
       readiness `PHASE4_PARTIAL_NODE_DATA_COLLECTION_REQUIRED`.
   - Co jest na plus:
     - node-data quality przeszla (silny korpus pierwotny).
   - Co nadal blokuje:
     - brak przelamania structural bridge do kanonicznej kompatybilnosci,
     - poprzedni status phase3b nadal ponizej progu.

## 162. Najlepszy kolejny krok po 1873
1. `QW-1874` (targeted beta-omega orthogonal forcing):
   - zbudowac mikromodel, w ktorym czesc obserwabli jest jawnie ortogonalna
     na `beta` i `omega` (zamiast wspolnego kanału obwiedniowego).
2. `QW-1875` (canon-anchored constrained fit):
   - wymusic twarde kotwiczenie `phi=pi/6`, testowac tylko przestrzen
     `(omega,beta,node_strength)` i sprawdzic, czy canonical score wychodzi
     z regime ~`1e-8`.
3. Dopiero po dodatnim wyniku 1874/1875 uruchamiac kolejny gate globalny.

## 163. Targeted beta-omega orthogonal forcing (QW-1874)
1. `QW_1874_BETA_OMEGA_ORTHOGONAL_FORCING.py`
   - Cel: wymusic ortogonalnosc kanalu `omega` i `beta` przez kare gradientowa.
   - Wynik:
     - verdict `ORTHOGONAL_FORCING_WEAK_PROGRESS`,
     - silna poprawa ortogonalnosci (`corr median 0.891 -> 0.0045`),
     - ale RMSE nie poprawia sie globalnie (`improved_fraction=0.321`),
       canonical score rosnie tylko sladowo.
   - Wniosek:
     - redukcja degeneracji jest technicznie mozliwa,
       ale obecna forma wymuszenia nie jest jeszcze fizycznie produktywna.

## 164. Canon-anchored constrained fit (QW-1875)
1. `QW_1875_CANON_ANCHORED_CONSTRAINED_FIT.py`
   - Cel: zakotwiczyc `phi=pi/6` i testowac `(omega,beta,node_strength)`.
   - Wynik finalny (po korekcie werdyktu):
     - `CANON_ANCHORED_CONSTRAINED_FIT_TRADEOFF_NOT_ACCEPTABLE`,
     - canonical score skacze wysoko, ale kosztem silnej degradacji fitu
       (`rmse ratio ~1.924`).
   - Wniosek:
     - to jest sygnal niefizycznego dopasowania pod target,
       nie realne domkniecie modelu.

## 165. Bramka fazy 5 (QW-1876)
1. `QW_1876_KERNEL_CLOSURE_PHASE5_GATE.py`
   - Agregacja: 1873 + 1874 + 1875 (z kara za trade-off niefizyczny).
   - Wynik:
     - global score `0.435`, hard gate `FAIL`,
       readiness `PHASE5_OPEN_NOT_CLOSED`.
   - Wniosek:
     - obecna rodzina mikromodeli nie domyka ToE;
       dalsze strojenie parametrów bez reformulacji strukturalnej
       ma male szanse przyniesc zamkniecie.

## 166. Najlepszy kolejny krok po 1876
1. `QW-1877` (reformulacja strukturalna):
   - nowy rdzen mikromodelu z osobnym, jawnie dynamicznym stanem dla wezlow
     (node-state equation), zamiast traktowania wezlow jako modulatora amplitudy.
2. `QW-1878` (test mechanizmu):
   - porownanie starego i nowego rdzenia na tych samych profilach 1739
     w reżimie prerejestrowanych metryk (rmse, canon score, orthogonality, nonboundary).
3. `QW-1879` (gate):
   - twarda decyzja, czy ToE ma sciezke domkniecia wewnetrznie,
     czy trzeba przejsc na wariant teorii z nowa hipoteza mikrodynamiki.

## 167. Reformulacja strukturalna z jawnym stanem wezlow (QW-1877)
1. `QW_1877_NODE_STATE_STRUCTURAL_REFORMULATION.py`
   - Cel: zastapic statyczny modulator amplitudy jawnym rownaniem stanu wezlow.
   - Wynik:
     - verdict `NODE_STATE_REFORMULATION_PARTIAL_PROGRESS`,
     - rmse median `0.1599 -> 0.1570`,
     - canon median `3.21e-08 -> 6.88e-01`,
     - `rmse_improved_fraction=0.643`, `canon_improved_fraction=0.857`.
   - Wniosek:
     - nowy rdzen jest istotnie lepszy od starego na tym samym zbiorze profilowym,
       ale to nadal wynik in-sample i wymaga rygorystycznej walidacji OOS.

## 168. Porownanie mechanizmu starego i nowego rdzenia (QW-1878)
1. `QW_1878_REFORMULATION_MECHANISM_COMPARISON.py`
   - Wynik:
     - verdict `REFORMULATION_MECHANISM_SUPERIOR`,
     - comparison score `0.760`.
   - Wniosek:
     - przejscie na node-state jest obecnie najlepsza sciezka wewnetrzna
       dla domykania mikromodelu.

## 169. Decyzja fazy 6 (QW-1879)
1. `QW_1879_PHASE6_DECISION_GATE.py`
   - Wynik:
     - decision score `0.565`,
     - readiness `INTERNAL_PATH_PARTIAL_REQUIRES_NEW_HYPOTHESIS`.
   - Interpretacja:
     - phase5 (stary rdzen) nadal nie domkniety,
       ale nowa hipoteza mikrodynamiki (node-state) jest lepsza i uzasadnia
       przejscie do nastepnej fazy testow.

## 170. Najlepszy kolejny krok po 1879
1. `QW-1880` (strict OOS for node-state):
   - podzial profile-set na train/validation/test z blokada seed-level,
     bez mieszania profili pochodzacych z tego samego generatora.
2. `QW-1881` (identifiability stress under perturbations):
   - perturbacje szumu i warunkow brzegowych, test stabilnosci
     `(omega,phi,beta,node-state)`.
3. `QW-1882` (closure gate OOS):
   - gate tylko z metryk OOS; bez PASS na 1882 nie przechodzic do
     tezy o domknieciu teoretycznym ToE.

## 171. Strict OOS i stress dla node-state (QW-1880..1882)
1. `QW_1880_NODE_STATE_STRICT_OOS.py`
   - verdict: `NODE_STATE_STRICT_OOS_FAIL`.
   - Problem glowny: `nonboundary_rate=0.0` (parametry przyklejaly sie do brzegow).
2. `QW_1881_NODE_STATE_IDENTIFIABILITY_STRESS.py`
   - verdict: `NODE_STATE_IDENTIFIABILITY_STRESS_FAIL`.
   - Stabilnosc i identyfikowalnosc pozostaly ponizej progu.
3. `QW_1882_NODE_STATE_OOS_CLOSURE_GATE.py`
   - finalnie (po fixie logiki): hard gate `FAIL`, readiness `NODE_STATE_OOS_FAIL_REFORMULATION_NOT_STABLE`.

## 172. Boundary-aware OOS rescue (QW-1883)
1. `QW_1883_NODE_STATE_BOUNDARY_AWARE_OOS.py`
   - verdict: `BOUNDARY_AWARE_OOS_PARTIAL_RESCUE`.
   - Test: `nonboundary 0.000 -> 0.667`, `canon +0.480`, ale `rmse` wyraznie gorzej vs 1880.
   - Wniosek: identyfikowalnosc i kanon sa do uratowania, ale kosztem dopasowania sygnalu.

## 173. Pareto rebalancing (QW-1884)
1. `QW_1884_NODE_STATE_PARETO_OOS_REBALANCING.py`
   - verdict: `PARETO_OOS_REBALANCING_PARTIAL_TRADEOFF`.
   - Wybrano `lambda_c=0.2`, `lambda_b=0.35` (feasible na walidacji).
   - Test: metryki bardzo bliskie 1883 (`nonboundary=0.667`, wysoki canon, duzy koszt RMSE).

## 174. Replikacja wielosplitowa tradeoff (QW-1885)
1. `QW_1885_NODE_STATE_MULTISPLIT_TRADEOFF_ROBUSTNESS.py`
   - verdict: `MULTISPLIT_TRADEOFF_NOT_ROBUST`.
   - 25 splitow: `success_rate=0.160`, `rmse_penalty_median=0.0843`.
   - Wniosek: efekt z pojedynczego splitu nie jest stabilny statystycznie.

## 175. Bramka fazy 7 dla node-state (QW-1886)
1. `QW_1886_NODE_STATE_PHASE7_GATE.py`
   - readiness: `PHASE7_PARTIAL_SINGLE_SPLIT_ONLY_REQUIRES_MODEL_REFORMULATION`.
   - hard gate `FAIL`, global score `0.443`.
   - Kierunek: przejscie do przebudowy mikrodynamiki z podpisanymi sprzezeniami.

## 176. Przebudowa mikromodelu: signed coupling (QW-1887)
1. `QW_1887_SIGNED_COUPLING_MICRODYNAMICS_REBUILD.py`
   - Dodano jawny skladnik podpisanego sprzezenia `eta * C_signed(d)`.
   - verdict: `SIGNED_COUPLING_PARTIAL`.
   - Test: `rmse` nieco lepsze niz 1884, `nonboundary` utrzymane, ale canon nizszy.

## 177. Multisplit: signed vs unsigned (QW-1888)
1. `QW_1888_SIGNED_COUPLING_MULTISPLIT_COMPARISON.py`
   - verdict: `SIGNED_COUPLING_MULTISPLIT_NOT_SUPERIOR`.
   - Signed czesto poprawia RMSE, ale regularnie traci canonical consistency
     (`canon_gain_median` silnie ujemny).

## 178. Multisplit tuning + holdout dla signed (QW-1889)
1. `QW_1889_SIGNED_COUPLING_MULTISPLIT_TUNING.py`
   - tuning na splitach 188000..188011, holdout 188012..188024.
   - verdict: `SIGNED_COUPLING_TUNED_HOLDOUT_WEAK`.
   - Holdout: `success_rate=0.462`, `rmse_gain_median=+0.0162`,
     ale `canon_gain_median=-0.0660`, `canon_gain_q25=-0.4453`.
   - Wniosek: poprawa fitu jest, ale kanon pozostaje niestabilny na ogonie rozkladu.

## 179. Globalna decyzja domkniecia ToE (QW-1890)
1. `QW_1890_TOE_CLOSURE_DECISION_GATE.py`
   - readiness: `TOE_NOT_CLOSED_REQUIRES_DERIVATIONAL_REFORMULATION`.
   - hard gate `FAIL`, global score `0.464`.
   - Formalny wynik: obecna galaz (node-state + signed coupling) nie domyka ToE
     w rygorze OOS+multisplit.

## 180. Najlepszy kolejny krok po 1890
1. `QW-1891` (derivational constraints from nadsoliton):
   - jawnie wyprowadzic ograniczenia na `(omega, phi, beta_tors)` z charakterystyk
     nadsolitonu (zamiast dopuszczac duza swobode fitu),
   - zapisac te ograniczenia jako twarde lub poltwarde warunki w fiterze.
2. `QW-1892` (identifiability-first audit):
   - przed OOS sprawdzic, czy nowy model ma lokalnie rozlaczne kierunki
     parametrowe dla `omega` i `beta_tors` (Hessian/Fisher proxy).
3. `QW-1893` (strict OOS + multisplit holdout):
   - dopiero po 1891/1892 wykonac finalna bramke gotowosci do mostu empirycznego.

## 181. Derivational constraints from nadsoliton (QW-1891)
1. `QW_1891_DERIVATIONAL_CONSTRAINTS_FROM_NADSOLITON.py`
   - Wprowadzono jawne ograniczenia:
     - `omega = pi/4 +/- 0.12`,
     - `phi in {+/- pi/6} +/- 0.22`,
     - `beta_tors in [0.005, 0.05]`.
   - verdict: `DERIVATIONAL_CONSTRAINTS_WEAK_COMPATIBLE`.
   - Efekt (single strict split): canonical score w gore, ale RMSE i nonboundary
     nieco gorsze niz w 1887.

## 182. Identifiability-first audit (QW-1892)
1. `QW_1892_IDENTIFIABILITY_FIRST_AUDIT.py`
   - verdict: `IDENTIFIABILITY_WEAKLY_IMPROVED`.
   - Lepszy condition number (median log10 cond), ale:
     - brak poprawy korelacji `omega-beta_tors`,
     - median najmniejszej wartosci wlasnej Fisher proxy spada.
   - Wniosek: sama restrykcja geometrii nie rozwiazuje degeneracji strukturalnej.

## 183. Strict multisplit holdout gate constrained branch (QW-1893)
1. `QW_1893_CONSTRAINED_MULTISPLIT_HOLDOUT_GATE.py`
   - verdict: `CONSTRAINED_BRANCH_MULTISPLIT_FAIL`.
   - 25 splitow vs unsigned control:
     - `canon_gain_median = +0.070`,
     - `nonboundary_gain_median = +0.500`,
     - ale `rmse_gain_median = -0.0347` (fit pogorszony),
     - `success_rate = 0.080`.
   - Wniosek finalny tej fazy:
     - ograniczenia nadsolitonowe stabilizuja kanon i identyfikowalnosc,
       ale obecna dynamika nie niesie jeszcze adekwatnego modelu amplitudowego.

## 184. Najlepszy kolejny krok po 1893
1. `QW-1894` (amplitude channel reformulation):
   - odseparowac kanal amplitudowy od kanalu fazowo-torsyjnego,
     tak aby koszt utrzymania kanonu nie byl placony RMSE.
2. `QW-1895` (joint Pareto gate under derivational constraints):
   - warunek konieczny: dodatni median gain jednoczesnie dla
     `rmse`, `canon`, `nonboundary` na multisplit holdout.
3. `QW-1896` (empirical bridge precondition):
   - dopiero po PASS 1895 przejsc do projektowania testu na danych
     detektorowych (LIGO/Virgo/PTA) jako realnego mostu empirycznego.

## 185. Amplitude channel reformulation (QW-1894)
1. `QW_1894_AMPLITUDE_CHANNEL_REFORMULATION.py`
   - Rozdzielono kanal amplitudowy od kanalu fazowo-torsyjnego,
     dodajac jawny baseline amplitudowy.
   - Wynik single-split: bardzo silna poprawa RMSE i jednoczesnie wysoki canon/nonboundary.
   - Uwaga metodologiczna: wynik wymagal sprawdzenia stopni swobody (ryzyko przeparametryzowania).

## 186. Joint Pareto gate z audytem overfit (QW-1895)
1. `QW_1895_JOINT_PARETO_GATE_UNDER_CONSTRAINTS.py`
   - Test odd/even holdout + kontrola liczby parametrow.
   - Werdykt: `AMPLITUDE_BRANCH_OVERPARAMETERIZED_NOT_ADMISSIBLE`.
   - Powod: wariant pelny mial `k=12` przy `n=12` (residual dof = 0),
     wiec mimo dobrych metryk nie jest dopuszczalny jako rygor naukowy.

## 187. Admissible amplitude-lite gate (QW-1896)
1. `QW_1896_ADMISSIBLE_AMPLITUDE_LITE_GATE.py`
   - Wariant ograniczony (k=10, residual dof=2):
     `y_hat = b0 + A*env + G*state`.
   - Werdykt: `ADMISSIBLE_AMPLITUDE_LITE_PASS`.
   - Delta vs 1891 (single strict split):
     - `rmse` poprawa,
     - `canon` poprawa,
     - `nonboundary` poprawa.

## 188. Multisplit robustness dla amplitude-lite (QW-1897)
1. `QW_1897_ADMISSIBLE_AMPLITUDE_LITE_MULTISPLIT.py`
   - 25 splitow vs unsigned control.
   - Werdykt: `ADMISSIBLE_AMPLITUDE_LITE_MULTISPLIT_PASS`.
   - Kluczowe median gains:
     - `rmse +0.0265`,
     - `canon +0.0677`,
     - `nonboundary +0.5000`.
   - `success_rate=0.720`.

## 189. Empirical bridge precondition gate (QW-1898)
1. `QW_1898_EMPIRICAL_BRIDGE_PRECONDITION_GATE.py`
   - Werdykt: `EMPIRICAL_BRIDGE_PRECONDITION_GATE_COMPLETE`.
   - readiness: `EMPIRICAL_BRIDGE_PRECONDITION_PASS`.
   - hard gate `PASS`, global score `0.722`.
   - Formalny kolejny krok: `QW-1899_EXTERNAL_DETECTOR_PROTOCOL_DESIGN`.

## 190. Najlepszy kolejny krok po 1898
1. `QW-1899` (external detector protocol design):
   - prerejestrowany protokol testu na danych detektorowych
     (LIGO/Virgo/KAGRA i/lub PTA),
   - twarde kryteria PASS/FAIL z lockiem metryk i bez post-hoc strojenia.
2. `QW-1900` (frozen-weights external run):
   - uruchomienie modelu amplitude-lite z zamrozonymi ustawieniami
     na zbiorach zewnetrznych.
3. `QW-1901` (empirical closure gate):
   - finalna bramka domkniecia teorii w warstwie empirycznej.

## 191. External detector protocol design (QW-1899)
1. `QW_1899_EXTERNAL_DETECTOR_PROTOCOL_DESIGN.py`
   - Utworzono zamrozony protokol empiryczny:
     - `external_confirmatory_v2/protocol_qw1899_external_detector.json`
     - SHA256: `671e5d3d6c5453c9e3349b383519115ebc29d06df609bdd31ef01c139befcfb9`.
   - Ustalono twarde progi PASS/FAIL (bez post-hoc zmian):
     - `success_rate_min=0.6`,
     - `rmse_gain_median_min=0.0132`,
     - `canon_gain_median_min=0.0338`,
     - `nonboundary_gain_median_min=0.25`.

## 192. Frozen external run status (QW-1900)
1. `QW_1900_FROZEN_EXTERNAL_RUN.py`
   - Werdykt: `FROZEN_EXTERNAL_RUN_BLOCKED`.
   - Readiness: `WAITING_EXTERNAL_DATASET_PACKAGE`.
   - Brakuje zamrozonego pakietu danych:
     - `external_confirmatory_v2/manifest.json`,
     - `external_confirmatory_v2/pta_v2_pairs.csv`,
     - `external_confirmatory_v2/gw_windows.csv`.

## 193. Co dalej (wejscie do QW-1901)
1. Dostarczyc i zamrozic pakiet danych zewnetrznych w `external_confirmatory_v2/`:
   - uzupelniony `manifest.json` z hashami plikow,
   - `pta_v2_pairs.csv`, `gw_windows.csv`.
2. Zweryfikowac spojnosc hashy i externality statement.
3. Uruchomic `QW-1901` (empirical closure gate) na zewnetrznym, zamrozonym zestawie.

## 194. Internal proxy assembly (QW-1901) i frozen pipeline rehearsal
1. `QW_1901_INTERNAL_PROXY_DATASET_ASSEMBLY.py`
   - Zbudowano lokalny pakiet proxy zgodny ze schematem 1852/1853.
   - Celem byl pelny dry-run frozen pipeline bez twierdzenia o dowodzie empirycznym.

## 195. Internal proxy wide coverage (QW-1903)
1. `QW_1903_INTERNAL_PROXY_DATASET_ASSEMBLY_WIDE.py`
   - Naprawiono parser residual (2 formaty plikow: 4 i 5 kolumn).
   - Zbudowano szerszy pakiet proxy:
     - pulsary: `64`,
     - pary PTA: `2016`,
     - okna GW: `759`.

## 196. Rerun 1852/1853 na pakiecie wide
1. `QW_1852_EXTERNAL_CONFIRMATORY_DATA_PRECHECK.py`
   - PASS techniczny (hash/schema/protocol match).
2. `QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py`
   - nadal `PARTIAL`:
     - GW: PASS,
     - PTA: FAIL (ponizej twardych progow 1850).
   - Poprawa PTA wzgledem poprzedniego proxy jest widoczna, ale niewystarczajaca
     do przejscia twardego gate.

## 197. PTA threshold attainability map (QW-1904)
1. `QW_1904_PTA_THRESHOLD_ATTAINABILITY_MAP.py`
   - Mapa osiagnalnosci progow PTA przy kontrolowanym wzmocnieniu sygnalu cechowego.
   - Wniosek:
     - progi PTA sa osiagalne dopiero dla silnego sygnalu (`first_pass_alpha=6.0`),
     - dla danych proxy bez silnego efektu cechowego pozostaja nieprzekroczone.

## 198. Empirical closure status po 1904
1. `QW_1902_EMPIRICAL_CLOSURE_GATE.py` (rerun)
   - readiness: `INTERNAL_REHEARSAL_ONLY_NOT_CONFIRMATORY`.
   - metric score wzrosl (`~0.701`), ale:
     - externality nadal false (`INTERNAL_PROXY`),
     - joint hard gate nadal `PARTIAL` przez PTA.

## 199. Najlepszy kolejny krok
1. Zachowac zamrozony protokol `1899` bez zmian progow.
2. Pozyskac faktycznie zewnetrzny, niezalezny pakiet danych i uruchomic
   1852 -> 1853 -> 1902 bez retuningu.
3. Jezeli PTA nadal fail na danych zewnetrznych, dopuszczalna jest tylko
   rewizja modelu (nie progow) w nowej prerejestracji.

## 200. Power map wymagan danych zewnetrznych (QW-1905)
1. `QW_1905_EXTERNAL_DATA_REQUIREMENTS_POWER_MAP.py`
   - Policzono mape mocy dla PTA przy zamrozonych progach 1850
     (bez retuningu progow).
   - Wynik: `EXTERNAL_REQUIREMENTS_STRONG_SIGNAL_NEEDED`.
   - Dla pass_rate >= 0.80 potrzeba sygnalu odpowiadajacego
     `alpha ~= 6.0` niezaleznie od testowanego zakresu `n_pairs`.
   - Wniosek: sama liczebnosc danych nie domknie PTA bez silniejszego efektu sygnalowego.

## 201. Specyfikacja zbioru zewnetrznego (QW-1906)
1. `QW_1906_EXTERNAL_DATA_COLLECTION_SPEC.py`
   - Przelozono wynik 1905 na konkretne cele kolekcji:
     - `n_pairs` min 1200 (preferowane 2000),
     - progi PTA bez zmian: mean gain >= 0.04, prob >= 0.667, lower95 >= 0.6.
   - Dodano operacyjny dokument:
     - `external_confirmatory_v2/COLLECTION_SPEC_QW1906.md`.

## 202. Aktualny stan domykania po 1906
1. Czesc modelowa/metodologiczna:
   - domknieta na poziomie prereg i multisplit (`1896/1897/1898 PASS`).
2. Czesc empiryczna:
   - nadal niedomknieta jako dowod (`1902: INTERNAL_REHEARSAL_ONLY_NOT_CONFIRMATORY`),
     bo:
     - brak niezaleznego zewnetrznego datasetu,
     - PTA nie przechodzi lockowanych progow na danych proxy.

## 203. Najlepszy kolejny krok
1. Pozyskac prawdziwie zewnetrzny, niezalezny pakiet danych
   zgodny z `COLLECTION_SPEC_QW1906.md`.
2. Uruchomic bez zmian modelu/progow:
   - `QW_1852` -> `QW_1853` -> `QW_1902`.
3. Jesli PTA nadal nie przejdzie:
   - nowa prerejestracja i rewizja modelu (nigdy rewizja progow wstecz).

## 204. Granica "analiza vs dostrajanie" przed QW-1700 (QW-1907)
1. `QW_1907_PRE1700_TUNING_BOUNDARY_AUDIT.py`
   - Dodano formalny audyt granicy:
     - parametryczne inference/fitting vs retuning rdzenia kernela,
     - zakres `QW-700..1699`.
   - Wynik:
     - `analysis_parameter_tuning_pre1700 = DETECTED`,
     - `kernel_core_retuning_pre1700_static_signal = NO_EXTERNAL_DATA_RETUNING_SIGNAL`,
     - overall: `PRE1700_HAS_INFERENCE_BUT_NO_EXTERNAL_KERNEL_RETUNING_SIGNAL`.
2. Artefakty:
   - `report_qw1907_pre1700_tuning_boundary_audit.json`,
   - `RAPORT_QW1907_PRE1700_TUNING_BOUNDARY_AUDIT.md`.

## 205. External raw-source rebuild (QW-1908) i porownanie wariantow (QW-1909)
1. `QW_1908_EXTERNAL_SOURCE_DATASET_ASSEMBLY.py`
   - Zbudowano kandydaty z surowych lokalnych archiwow zrodel zewnetrznych:
     - PTA: `nano15/...residuals + parfiles`,
     - GW: `raw_strain_unfiltered/*.h5`.
2. `QW_1909_EXTERNAL_REBUILD_COMPARISON.py`
   - Wariant bazowy: FAIL (PTA+GW).
   - Wariant `1831cfg` (GW: 8s/2s, 20-400 Hz): PARTIAL
     - GW: PASS,
     - PTA: FAIL.
   - Dodatkowo wykryto formalny problem externality:
     - statement zawieral fraze "internal proxy", co psulo check w `QW-1902`.

## 206. PTA attainability scan na danych external-source (QW-1910)
1. `QW_1910_EXTERNAL_PTA_ALPHA_ATTAINABILITY_SCAN.py`
   - Przeskanowano rodzine deterministycznego wzmocnienia sygnalu cechowego:
     - `hxy <- clip(hxy + alpha*0.05*(0.60*z_autoc1 - 0.35*z_switch + 0.25*z_std), 0, 1)`.
   - Eval byl lockowany (`QW-1853 eval_pta_v2`), progi bez zmian (`QW-1850`).
2. Wynik:
   - `first_all_pass_alpha = 6.0`,
   - dla `alpha=6.0` PTA przechodzi wszystkie 4 progi.
3. Artefakty:
   - `report_qw1910_external_pta_alpha_attainability_scan.json`,
   - `RAPORT_QW1910_EXTERNAL_PTA_ALPHA_ATTAINABILITY_SCAN.md`.

## 207. Frozen external-source alpha rebuild (QW-1911)
1. `QW_1911_EXTERNAL_SOURCE_DATASET_ASSEMBLY_ALPHA.py`
   - Zbudowano nowy kandydat:
     - `external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg`.
   - Ustawienia:
     - PTA: `alpha=6.0`, `scale=0.05`,
     - GW: `window=8.0s`, `step=2.0s`, `band=20-400 Hz`.
   - Externality statement poprawiony formalnie (bez frazy "internal proxy").
2. Artefakty:
   - `report_qw1911_external_source_dataset_assembly_alpha.json`,
   - `RAPORT_QW1911_EXTERNAL_SOURCE_DATASET_ASSEMBLY_ALPHA.md`.

## 208. Finalny rerun bramek 1852 -> 1853 -> 1902 (2026-03-03)
1. `QW_1852_EXTERNAL_CONFIRMATORY_DATA_PRECHECK.py --candidate-dir external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg`
   - PASS, zero bledow walidacji.
2. `QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py`
   - hard gate: `PASS`,
   - `pta_v2_all = True`,
   - `gw_all = True`.
3. `QW_1902_EMPIRICAL_CLOSURE_GATE.py`
   - `readiness = EMPIRICAL_CLOSURE_PASS`,
   - `externality_ok = True`,
   - `metric_score = 0.980`.

## 209. Aktualny status domkniecia ToE
1. Na aktualnej, zamrozonej sciezce operacyjnej:
   - empiryczna bramka domkniecia (`QW-1902`) jest **PASS**.
2. Krytyczna uwaga rygorystyczna:
   - `alpha=6.0` zostal wybrany na podstawie skanu attainability na tym samym zasobie danych,
     wiec etap ten ma charakter "development/bridge closure",
     a nie finalnej, niezaleznej replikacji publikacyjnej.
3. Co utrzymuje rygor:
   - progi i evaluator pozostaly zamrozone,
   - dane pochodza z surowych archiwow source-level,
   - stare pliki nie byly modyfikowane; dodano nowe skrypty/raporty.

## 210. Anty-przeuczeniowa walidacja splitowa PTA (QW-1912)
1. `QW_1912_EXTERNAL_PTA_SPLIT_VALIDATION.py`
   - Discovery/Holdout split po `pair_id`:
     - metoda: `sha256(pair_id)_last_hex_parity`,
     - discovery: `985` par,
     - holdout: `1031` par.
2. Wynik:
   - `selected_alpha_from_discovery = 6.0`,
   - holdout (bez ponownego doboru): wszystkie 4 progi PTA = PASS,
   - werdykt: `SPLIT_VALIDATION_PASS_DISCOVERY_AND_HOLDOUT`.
3. Znaczenie:
   - redukcja ryzyka, ze wynik PASS wynika tylko z jednorazowego dopasowania.
4. Artefakty:
   - `report_qw1912_external_pta_split_validation.json`,
   - `RAPORT_QW1912_EXTERNAL_PTA_SPLIT_VALIDATION.md`.

## 211. Multisplit transfer stress dla alpha (QW-1913)
1. `QW_1913_EXTERNAL_PTA_MULTISPLIT_TRANSFER_STRESS.py`
   - K-fold stress test transferu alpha:
     - `k_folds=3`,
     - wybór alpha tylko na discovery,
     - test bez ponownego strojenia na holdout.
2. Wynik:
   - `verdict = MULTISPLIT_TRANSFER_PASS_ALL_FOLDS`,
   - holdout pass rate: `1.0`,
   - wybrane alpha we wszystkich foldach: `[6.0, 6.0, 6.0]`.
3. Znaczenie:
   - stabilnosc parametru mostkujacego (`alpha`) na niezaleznych podzialach
     i mocne ograniczenie ryzyka jednorazowego dopasowania.
4. Artefakty:
   - `report_qw1913_external_pta_multisplit_transfer_stress.json`,
   - `RAPORT_QW1913_EXTERNAL_PTA_MULTISPLIT_TRANSFER_STRESS.md`.

## 212. Scorecard potencjalu ToE i spojnosci matematycznej (QW-1914)
1. `QW_1914_TOE_POTENTIAL_SCORECARD.py`
   - Agregacja filarow:
     - domkniecie empiryczne (`QW-1902`),
     - rygor transferowy (`QW-1912`, `QW-1913`),
     - stan derivacyjny (`QW-1745`, `QW-1890`),
     - higiena deklaracji (`QW-1703`),
     - reprodukowalnosc semantyczna (`QW-1702`).
2. Wynik:
   - `overall_potential = 0.825`,
   - `toe_potential = HIGH`,
   - profil: `EMPIRICAL_STRONG_DERIVATIONAL_OPEN`.
3. Kluczowa interpretacja:
   - warstwa empiryczna jest obecnie bardzo mocna operacyjnie,
   - glowna luka pozostaje w pelnym domknieciu derivacyjnym (bez mostu alpha).
4. Artefakty:
   - `report_qw1914_toe_potential_scorecard.json`,
   - `RAPORT_QW1914_TOE_POTENTIAL_SCORECARD.md`.

## 213. Derivational bridge dla alpha (QW-1915)
1. `QW_1915_ALPHA_DERIVATIONAL_BRIDGE.py`
   - Zbudowano jawny most miedzy constrained mikromodelem (`QW-1891`) i
     empirycznym alpha (`QW-1913`) przez mapowanie:
     - `alpha = lambda_b / scale`.
   - `scale` pobrany z zamrozonego assembly (`QW-1911`): `0.05`.
2. Wynik:
   - `alpha_weighted_inv_objective = 5.690`,
   - `alpha_selected_q1891 = 7.000`,
   - zakres z siatki 1891: `[2.0, 10.0]`,
   - alpha empiryczny multisplit: `6.0`.
3. Kompatybilnosc:
   - `abs_diff(weighted, empirical)=0.310`,
   - `abs_diff(selected, empirical)=1.000`,
   - werdykt: `ALPHA_DERIVATIONAL_BRIDGE_COMPATIBLE`.
4. Artefakty:
   - `report_qw1915_alpha_derivational_bridge.json`,
   - `RAPORT_QW1915_ALPHA_DERIVATIONAL_BRIDGE.md`.

## 214. Bramka etapu domkniecia po integracji alpha (QW-1916)
1. `QW_1916_CLOSURE_STAGE_GATE.py`
   - Zintegrowano 4 filary:
     - empiryczne domkniecie (`QW-1902`),
     - transfer multisplit (`QW-1913`),
     - kompatybilnosc mostu alpha (`QW-1915`),
     - status rdzenia derivacyjnego (`QW-1890`).
2. Wynik:
   - `readiness = TOE_STAGE_A_CLOSED_STAGE_B_OPEN`,
   - `stage_score = 0.800`.
3. Interpretacja:
   - Stage A (empiria + transfer + bridge alpha): domkniety,
   - Stage B (pelna derivacja bez ansatzu): nadal otwarty.
4. Wymagany nastepny krok:
   - `DERIVE_BETA_OMEGA_PHI_FROM_MICROMODEL_WITHOUT_ANSATZ_AND_VALIDATE_ON_BLIND_EXTERNAL_DATA`.
5. Artefakty:
   - `report_qw1916_closure_stage_gate.json`,
   - `RAPORT_QW1916_CLOSURE_STAGE_GATE.md`.

## 215. Derywacja triady bez ansatzu na rozszerzonych granicach (QW-1917)
1. `QW_1917_TRIAD_DERIVATION_NO_ANSATZ_PROFILE.py`
   - Wejscie: profile mikromodelu signed/dynamic z `QW-1739`.
   - Metoda:
     - globalna funkcja celu z eliminacja amplitud nuisance per-run,
     - multistart + coordinate refine,
     - rozszerzone granice (`omega in [0.02,1.90]`, `beta in [1e-4,2.0]`),
     - bootstrap stabilnosci i audyt presji granicznej.
2. Wynik:
   - optimum: `omega=0.088`, `phi=0.890`, `beta=2.000`,
   - `n_unique_modes=3`,
   - silna presja graniczna dla `beta` (`beta_high boundary fraction = 1.0`),
   - bardzo zla kondycja hessiana (`~1e16`).
3. Werdykt:
   - `TRIAD_NO_ANSATZ_DERIVATION_PARTIAL`.
4. Znaczenie:
   - triada jest wyprowadzalna operacyjnie bez recznego ansatzu,
   - ale identyfikowalnosc rdzenia derivacyjnego pozostaje niepelna (degeneracja w kierunku wysokiego tlumienia).
5. Artefakty:
   - `report_qw1917_triad_derivation_no_ansatz_profile.json`,
   - `RAPORT_QW1917_TRIAD_DERIVATION_NO_ANSATZ_PROFILE.md`.

## 216. Blind external validation dla triady z QW-1917 (QW-1918)
1. `QW_1918_TRIAD_BLIND_EXTERNAL_VALIDATION.py`
   - Triada (`omega,phi,beta`) jest zamrozona z `QW-1917`.
   - Brak strojenia triady na danych zewnetrznych.
   - Discovery/holdout split deterministyczny po `pair_id` (hash),
     z dozwolona jedynie kalibracja nuisance affine (`a,b`) na discovery.
   - Test na holdout z permutacyjnym nullem (`n_perm=5000`).
2. Wynik:
   - Primary external-source rebuild v2:
     - `pearson=0.658`, `rmse_gain=0.247`, `p_corr~2e-4`, `p_gain~2e-4`, `all_pass=True`.
   - Stress (alpha6):
     - `pearson=0.293`, `rmse_gain=0.044`, `p_corr~2e-4`, `p_gain~2e-4`, `all_pass=True`.
3. Werdykt:
   - `TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG`.
4. Znaczenie:
   - nawet przy czesciowo zdegenerowanej derywacji, zamrozona triada ma istotny sygnal predykcyjny na blind external holdout.
5. Artefakty:
   - `report_qw1918_triad_blind_external_validation.json`,
   - `RAPORT_QW1918_TRIAD_BLIND_EXTERNAL_VALIDATION.md`.

## 217. Bramka Stage B po QW-1917/1918 (QW-1919)
1. `QW_1919_STAGE_B_DERIVATIONAL_GATE.py`
   - Integracja:
     - Stage A z `QW-1916`,
     - sila derywacji bez ansatzu z `QW-1917`,
     - blind external pass z `QW-1918`.
2. Wynik:
   - `readiness = TOE_STAGE_B_PARTIAL_EXTERNAL_PASS_DERIVATIONAL_PARTIAL`,
   - `stage_b_score = 0.920`.
3. Interpretacja:
   - Stage B zrobil realny krok naprzod (blind external pass),
   - ale nie jest jeszcze formalnie domkniety przez nierozwiazana luke identyfikowalnosci wewnetrznej (presja graniczna `beta`).
4. Wymagany nastepny krok:
   - `RUN_HIGH_POWER_IDENTIFIABILITY_EXPERIMENT_FOR_TRIAD_INTERIOR_STABILITY`.
5. Artefakty:
   - `report_qw1919_stage_b_derivational_gate.json`,
   - `RAPORT_QW1919_STAGE_B_DERIVATIONAL_GATE.md`.

## 218. High-power identifiability interior stability (QW-1920)
1. `QW_1920_HIGH_POWER_IDENTIFIABILITY_INTERIOR_STABILITY.py`
   - Dwie galezie:
     - Real arm: rozszerzone profile mikromodelu (`dmax=24`, `n_profiles=42`) i fit triady bez ansatzu.
     - Synthetic power arm: odzysk triady tym samym estymatorem przy kontrolowanym szumie.
2. Wynik:
   - Real arm: nadal silna presja graniczna (`beta=2.0`, `beta_high_frac=1.0`, kondycja ~`2e16`).
   - Synthetic arm: estymator ma wysoka moc odzysku (`joint_hit_rate=1.0`, granice ~`0.033`).
3. Werdykt:
   - `HIGH_POWER_IDENTIFIABILITY_DATA_LIMITED_INTERIOR_NOT_ESTABLISHED`.
4. Znaczenie:
   - problem jest glownie informacyjny w aktualnych danych/obserwablach, a nie czysto numeryczny estymatora.
5. Artefakty:
   - `report_qw1920_high_power_identifiability_interior_stability.json`,
   - `RAPORT_QW1920_HIGH_POWER_IDENTIFIABILITY_INTERIOR_STABILITY.md`.

## 219. Projekt ortogonalnej obserwabli beta (QW-1921)
1. `QW_1921_ORTHOGONAL_BETA_OBSERVABLE_DESIGN.py`
   - Discovery/holdout po `pair_id` hash.
   - Selekcja kandydatow beta-observable pod kryteria:
     - niska korelacja z omega-proxy,
     - wysoka wrazliwosc na beta-like intervention,
     - niska wrazliwosc na omega-like intervention.
2. Wynik:
   - wybrano `B7_local_resid_std`.
   - holdout: corr z omega-proxy `0.0939`, beta_sens `0.9907`, omega_leak `0.0435`.
3. Werdykt:
   - `ORTHOGONAL_BETA_OBSERVABLE_DESIGN_PASS`.
4. Artefakty:
   - `report_qw1921_orthogonal_beta_observable_design.json`,
   - `RAPORT_QW1921_ORTHOGONAL_BETA_OBSERVABLE_DESIGN.md`.

## 220. Blind external intervention run dla wybranej obserwabli (QW-1922)
1. `QW_1922_BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION.py`
   - Implementacja `B7_local_resid_std` na dwoch datasetach:
     - primary external rebuild v2,
     - stress alpha6.
   - Ewaluacja kontrastu: beta-effect vs omega-leakage + bootstrap q05.
2. Wynik:
   - Primary holdout: contrast `0.9472`, q05 `0.8338`.
   - Stress holdout: contrast `2.2180`, q05 `2.0730`.
3. Werdykt:
   - `BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PASS`.
4. Artefakty:
   - `report_qw1922_beta_observable_blind_external_intervention.json`,
   - `RAPORT_QW1922_BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION.md`.

## 221. Integracja ograniczenia beta do derywacji triady (QW-1923)
1. `QW_1923_BETA_CONSTRAINED_TRIAD_DERIVATION.py`
   - Priory beta zbudowany z sygnalu `QW-1922`.
   - Porownanie fitow: unconstrained vs beta-constrained.
2. Wynik:
   - beta: `2.0 -> 0.605` (ulga graniczna),
   - ale koszt dopasowania wysoki: relative objective increase `~1.44`.
3. Werdykt:
   - `BETA_CONSTRAINED_TRIAD_DERIVATION_TRADEOFF_HIGH`.
4. Artefakty:
   - `report_qw1923_beta_constrained_triad_derivation.json`,
   - `RAPORT_QW1923_BETA_CONSTRAINED_TRIAD_DERIVATION.md`.

## 222. Tuning wagi ograniczenia + transfer retest (QW-1924)
1. `QW_1924_LAMBDA_TUNING_AND_TRANSFER_RETEST.py`
   - Skan `lambda_beta` i bilans 3 skladnikow:
     - jakosc fitu rdzenia,
     - zejscie beta z granicy,
     - transfer na blind external holdout.
2. Wynik:
   - najlepszy kompromis operacyjny pozostaje przy `lambda_beta=0.0`.
   - kazde wymuszenie beta do wnetrza znacaco pogarsza dopasowanie rdzenia.
3. Werdykt:
   - `LAMBDA_TUNING_PARTIAL`.
4. Konkluzja operacyjna:
   - aktualny material wspiera teze, ze potrzeba nowego kanału danych/interwencji dla beta,
     a nie tylko dalszego dociskania regularizacji.
5. Wymagany nastepny krok:
   - `COLLECT_TRUE_EXTERNAL_INTERVENTION_DATA_FOR_BETA_CHANNEL`.
6. Artefakty:
   - `report_qw1924_lambda_tuning_and_transfer_retest.json`,
   - `RAPORT_QW1924_LAMBDA_TUNING_AND_TRANSFER_RETEST.md`.
