# Jądro FIN jako operator propagacji informacji — niezależna analiza matematyczna

**Przedmiot:** repozytorium `github.com/hyconiek/Fractal-Nadsoliton-Theory`, stan na HEAD `db9a15d3` (2026-07-16)
**Pytanie badawcze:** *Czy matematyczna struktura jądra FIN może być interpretowana jako adaptacyjny operator propagacji informacji, którego naturalnym rozwinięciem byłaby samoucząca się sieć dynamiczna (np. hebbowska), a nie klasyczna teoria pola?*
**Metoda:** niezależna weryfikacja — żadne twierdzenie z audytów repo nie jest przyjmowane jako fakt; wszystkie wyniki spektralne i numeryczne poniżej policzone od zera z definicji jąder wydobytych z kodu źródłowego. Perspektywa „Teorii Wszystkiego" celowo pominięta.

---

## 0. Analiza ostatnich commitów (stan na 2026-07-16)

Ostatnie ~15 commitów (2026-06-26 → 2026-07-16) nie dotykają matematyki jądra. To wyłącznie **księgowość meta-poziomu** programu „strict-core closure":

| Commit | Data | Treść | Istota |
|---|---|---|---|
| `db9a15d3` | 07-16 | **Recenzja adwersarialna + roadmapa FIN** (`Recenzja_adwersarialna_i_roadmapa_FIN.md`) | dokument stwierdza m.in.: jądro to ansatz na losowym grafie geometrycznym; dynamika `dψ/dt = iKψ`; zero wielkości fizycznych wyprowadzonych ze „ścisłego rdzenia" |
| `8ab77781` (merge PR#302), `f016ad71` | 07-16 | **P3168–P3170**: certyfikat braku jednostki strict, audyt źródła ternarnego, twierdzenie o hitting-secie S₊ | `accepted_S_plus_sources: 0` — żaden kandydat na „źródło" nie przechodzi bramek |
| `64da1401` | 07-16 | `SUMMARY_GROK.md` | kolejne streszczenie AI |
| `3e51b9dd`, `2a2a68aa`, `c243e382` | 07-13/14 | wyczerpanie źródeł monomialnych S₊; audyt dopasowania prędkości światła na brzegu | zamknięcia dróg, brak konstrukcji |
| `466f9516`, `03722022`, `0f7e57d0`, `007fc983` | 07-13 | P3158 (rekonsyliacja zależności źródła jednostki), P3150 (audyt promienia hiperładunku) | obstrukcje |
| `11de33a6`, `a8566ff5` | 07-13 | **P3145**: „strict kernel SM GR reverse layout" | próba odwrócenia łańcucha jądro→SM/GR — niedomknięta |
| `cd539462`, `98a9e0b2` | 07-13 | P3141: audyt potencjału selektora aksjomatów | obstrukcja selektora |
| `4d606d9b`, `27cfeb2e` | 07-12 | **P3137**: audyt źródła ramy Fouriera (`fourier_frame_source_law`) | **wynik: `candidates passing all gates: 0`** — trybiki Fouriera Z₁₂ istnieją (strona „odbiorcza"), ale żadne źródło nie wybiera prymitywnego charakteru bez importowania |
| `0663b5d1` | 07-12 | P3132/S2082: guardrail przekroju lokalnego | meta |

**Diagnoza trajektorii:** projekt od miesięcy nie rozwija obiektu matematycznego (jądra), tylko dokumentuje braki własnej aksjomatyki („źródła", selektory, jednostki). Jądro strict jest **zamrożone** od fazy QW-2118 i występuje w nowych commitach wyłącznie jako stała referencja. Jedyny nowy dokument dotykający istoty pytania badawczego to recenzja adwersarialna z `db9a15d3`, która (§4.1) *przypuszcza*, że K(d) „wygląda jak funkcja Greena równania Helmholtza z tłumieniem" — bez dowodu. Własny, starszy audyt repozytorium (QW-2139) stwierdza natomiast wprost: *„duże residua radialne Laplace/Helmholtz (>1) odrzucają klasyczny lokalny Green claim"*. Poniżej rozstrzygam tę kwestię niezależnie i ostrzzej niż oba dokumenty.

---

## CZĘŚĆ I — IDENTYFIKACJA OPERATORA

### 1.1. Co faktycznie jest w kodzie (definicje wydobyte ze źródeł)

Dwie zamrożone postaci jądra (potwierdzone w `RELEASE_6_2_SHANNON_VOID_ASYMMETRY_EN_PL.md`, `QW_2133`, `report_qw2118`):

```
K_legacy_ont(d)  = α_geo · cos(ω d + φ) / (1 + β_tors · d)
                   α_geo = 4 ln 2 = 2.772589,  ω = π/4 ≈ 0.7854,  φ = π/6 ≈ 0.5236,  β_tors = 0.01

K_strict_gate(d) = cos(ω d + φ) / (1 + β · d^η)
                   ω = 0.18575,  φ = 0.1625,  β = 1.0,  η = 1.8
```

Dziedzina operacyjna: cykl Z₁₂ (12 „oktaw"), odległość cykliczna `d = min(|i−j|, 12−|i−j|) ∈ {1,…,6}`, zerowa diagonala; dynamika `dψ/dt = iKψ` (np. `QW-593`). Jądro występuje też na losowych chmurach punktów w ℝ³ (euclidesowe macierze losowe).

### 1.2. Czy to funkcja Greena Helmholtza? — wynik niezależnej analizy

**Nie — w żadnym z trzech sensownych reżimów. Ale odpowiedź jest ciekawsza niż samo „nie".**

**(a) Jako funkcja radialna na ℝ³ (interpretacja continuum).**
Funkcja Greena wychodząca Helmholtza 3D: `G⁺(r) = e^{iωr}/(4πr)` — zespolona, z fasz konkretną. Jądro legacy jest rzeczywiste i asymptotycznie (dla r ≫ a = 1/β = 100):

```
K_legacy(r) ~ (α_geo/β)·cos(ωr+φ)/r = (4πα_geo/β)·[cos φ·Re G⁺(r) − sin φ·Im G⁺(r)]
            = (4πα_geo/β)·[0.866·Re G⁺(r) − 0.500·Im G⁺(r)]
```

czyli **fala stojąca**: mieszanina części rzeczywistej rezolwenty (86,6%) i **miary spektralnej** operatora −Δ na powłoce q = ω (50%). Wychodząca funkcja Greena to kombinacja `cos + i·sin` z równymi wagami i fazą i; jądro FIN nie jest i nie może być (jest rzeczywiste) wychodzącą funkcją Greena. Kluczowe jednak: **w całym zakresie operacyjnym (d ≤ 11) tłumienie legacy wynosi < 9%** (obwiednia 2,745 → 2,498). Jądro legacy jest tam **prawie czystą falą stojącą o stałej amplitudzie** — czyli obiektem *on-shell* (funkcją własną / miarą spektralną na powłoce), a nie rezolwentą (obiektem *off-shell*). Obwiednia typu 1/r pojawia się dopiero przy r ~ 100, poza dziedziną użycia.

**(b) Na cyklu Z₁₂ (rzeczywista dziedzina).** Tu wszystko jest dokładne (brak przybliżeń). Oba jądra są cyrkulantami ⇒ diagonalizuje je DFT; symbole `λ_m = 2Σ_{d=1..5} K(d)cos(2πmd/12) + K(6)(−1)^m` w funkcji wartości własnych laplasjanu cyklu `μ_m = 2−2cos(2πm/12)`:

| m | μ_m | λ_m (legacy)/α_geo | λ_m (strict) |
|---|---|---|---|
| 0 | 0.000 | **−4.030** | **+1.660** |
| 1 | 0.268 | +0.742 | +0.906 |
| 2 | 1.000 | **+3.680** | +0.083 |
| 3 | 2.000 | −1.157 | −0.301 |
| 4 | 3.000 | −0.091 | −0.539 |
| 5 | 3.732 | −1.000 | −0.638 |
| 6 | 4.000 | −0.318 | −0.682 |

Rezolwenta `(zI − L)⁻¹` dowolnego operatora o wspólnej bazie Fouriera ma symbol **monotoniczny w μ_m z jednym biegunem** (jedna zmiana znaku, |λ| rośnie w stronę bieguna). Tymczasem:
- **legacy**: wzorzec znaków `− + + − − − −` — dwie zmiany znaku ⇒ **nie jest rezolwentą żadnego operatora jednobiegunowego**; jest **filtrem pasmowym** z pasmem wokół m = 1–2, czyli dokładnie wokół powłoki k ≈ ω = π/4 ⇒ **projektor spektralny na powłokę ω** (obiekt on-shell) — potwierdza wniosek (a);
- **strict**: wzorzec `+ + + − − − −` z |λ| *malejącym* w stronę przekroczenia zera ⇒ też nie rezolventa. Natomiast dopasowanie `λ_m ≈ c·(z − μ_m)·e^{−t μ_m}` daje **z = 1,179, t = 0,456, korelację 0,993** (maks. błąd wzgl. 8,8%). To znaczy:

> **K_strict jest spektralnie prawie dokładnie zgładzonym (pasmowo ograniczonym) operatorem Helmholtza `c·(z − L)·e^{−tL}` — czyli jest bliższy SAMEMU OPERATOROWI niż jego funkcji Greena.** Funkcją Greena byłoby K⁻¹, nie K.

To odwraca tezę z audytów: repo twierdzi „jądro ≈ funkcja Greena", a jest „jądro ≈ (zgładzony) operator; rezolwenta jest odwrotnością". Zresztą konsekwentnie: w dynamice `ψ̇ = iKψ` jądro gra rolę **generatora (hamiltonianu)**, nie propagatora — propagatorem jest `e^{iKt}`.

**(c) Jako jądro 1D.** Funkcja Greena Helmholtza 1D ma stały moduł (`e^{iω|x|}/2iω`). Strict zanika jak d^−1,8 ⇒ nie jest funkcją Greena 1D. Ciekawostka: obwiednia dalekiego pola Helmholtza w wymiarze d to r^{−(d−1)/2}; η = 1,8 odpowiadałoby `d_eff = 4,6` — nie ma wymiaru, w którym strict jest funkcją Greena (d = 4,6 ∉ ℕ).

### 1.3. Odpowiedzi na pytania szczegółowe Części I

- **Czy to funkcja Greena Helmholtza?** Nie. **Tylko asymptotycznie i tylko częściowo**: legacy ma asymptotycznie (r ≫ 100) obwiednię 3D Helmholtza, lecz jako falę stojącą (mieszanina Re G i Im G), nie wychodzącą funkcję Greena.
- **Dokładniejsza interpretacja:** legacy = **miara spektralna / rzut spektralny na powłokę q = ω** operatora typu −Δ (obiekt on-shell; w fizyce: „gęstość stanów na powłoce masowej", nie propagator). Strict = **operator różniczkowy eliptyczny pierwszego rzędu w L z wycięciem cieplnym w UV** (obiekt off-shell), sam będący „Helmholtzianem", którego funkcja Greena byłaby K⁻¹.
- **Bliższy inny operator?** Tak: strict leży bliżej rodziny `{(z−Δ)·e^{−tΔ}}` (operatory screened-Poisson/Helmholtz z regularyzacją) niż rodziny rezolwent. Legacy leży najbliżej **jądra miary spektralnej** `δ(−Δ − ω²)` (w 3D: `sin(ωr)/r` po zdjęciu fazy φ).
- **Symbol Fouriera:** policzony dokładnie powyżej (tabela). Na ℝ³: transformata Hankela ma osobliwość przy q = ω typu biegun PV (składnik cos) + skok/log (składnik sin) — struktura powłokowa, z UV zmienionym przez regularyzację `1+βd`.
- **Przestrzeń własna:** na cyklu — tryby Fouriera (dekalogiczne dublety m, 12−m ⇒ fale stojące cos/sin). Na losowej chmurze punktów — spektrum euklidesowej macierzy losowej (znana literatura: Mezzadri; Bogomolny–Bohigas–Schmit; Goetschy–Skipetrov).
- **Efektywna ranga** (iloraz partycypacji |λ|): **legacy 6,36 / 12; strict 8,21 / 12**. Ale koncentracja: legacy ma 43,6% masy widma w 2 modach (dublet m=2) — efektywnie „prawie ranga 2" w górnej części widma; strict jest widmowo znacznie bardziej rozmaity.
- **Najważniejsze własności spektralne:** (i) **żadne jądro nie jest dodatnio określone** (min λ: legacy −11,17; strict −0,68) ⇒ nie są jądrami Mercera ⇒ brak RKHS, brak interpretacji kowariancji GP; (ii) legacy — silna koncentracja na powłoce ω; (iii) strict — **elementowo nieujemne** (0% ujemnych sprzężeń) ⇒ korzeń Perrona–Frobeniusa λ₀ = 1,660 z dodatnim wektorem własnym: legalna macierz sąsiedztwa grafu ważonego mimo nieokreśloności spektralnej; (iv) sumy wierszowe strict stałe (cyrkulant) ⇒ tryb zerowy to wektor jednorodny.

---

## CZĘŚĆ II — K_legacy_ont KONTRA K_strict_gate

W repo występują jeszcze inne warianty historyczne (4-składnikowe `K_total = K_geo·K_res·(1+0,2·K_tors)·K_topo` z okresu QW-671, z tłumieniem wykładniczym zastąpionym hiperbolicznym „przez sumowanie ścieżek fraktalnych"), ale wszystkie są przodkami legacy i nie wnoszą nowej struktury operatorowej. Analizuję dwa kanoniczne.

### 2.1. Czym różnią się matematycznie

| Cecha | K_legacy_ont | K_strict_gate |
|---|---|---|
| Częstość ω | π/4 ≈ 0,785 (szybka oscylacja; okres 8 w jednostkach d) | 0,18575 (wolna; okres ≈ 33,8 ≫ średnica grafu) |
| Faza φ | π/6 ≈ 0,524 | 0,1625 (prawie czysty cosinus) |
| Tłumienie | 1/(1+0,01d) — prawie brak w dziedzinie (−9% na d=1…11) | 1/(1+d^1,8) — silne (÷38 na d=1…11) |
| Amplituda | α_geo = 4 ln 2 (skala) | 1 (bez prefaktora) |
| Znak sprzężeń | 72,7% ujemnych | 100% nieujemnych |
| Lokalność (udział masy wiersza przy d=1 / d≤2) | 8,5% / 24,9% — **silnie nielokalne** | 56,6% / 79,7% — **lokalne** |
| Symbol λ(μ) | pasmowy, on-shell (2 zmiany znaku) | filtrowany operator (z−μ)e^{−tμ}, korelacja dopasowania 0,993 |
| PSD / Mercer | nie (λ_min = −11,17) | nie (λ_min = −0,68), ale elementowo ≥ 0 (Perron–Frobenius) |
| Efektywna ranga | 6,36 (koncentracja: 43,6% w dublecie m=2) | 8,21 |

### 2.2. Odpowiedzi na pytania Części II

- **Które jest bliższe funkcji Greena?** Legacy — i to wyłącznie w sensie asymptotycznej obwiedni ℝ³ (1/r, r ≫ 100). To jednak „bliskość" kosmetyczna: obiekt on-shell. Strict nie jest blisko żadnej funkcji Greena — za to **jego odwrotność K_strict⁻¹ jest z definicji funkcją Greena operatora K_strict**, a ponieważ K_strict ≈ c(z−L)e^{−tL}, więc K_strict⁻¹ ≈ e^{tL}(z−L)⁻¹/c — rezolwenta Helmholtza z cieplną korektą UV. **To jest ścisła, poprawna wersja intuicji repo: nie „K jest funkcją Greena", lecz „K_strict jest (zgładzonym) operatorem, którego rezolwenta jest funkcją Greena".**
- **Które lepiej zachowuje własności operatora Helmholtza?** Strict — jego symbol jest pierwszego rzędu w μ z jednym zerem (z = 1,18 między μ₂ i μ₃), jak (z − Δ). Legacy ma symbol pasmowy, nemonotoniczny, dwuznakowy — struktura obca operatorom różniczkowym, naturalna dla projektorów spektralnych.
- **Które lepiej zachowuje lokalność?** Strict, zdecydowanie (80% masy w sąsiedztwie d ≤ 2). Legacy jest kwintesencjalnie nielokalne: prawie wszystkie węzły sprzęgają się porównywalnie silnie.
- **Które lepiej zachowuje interpretację propagatora?** Strict: elementowa nieujemność daje legalną interpretację wag przejść / adjacency (a po normalizacji wierszowej — operatora przejścia ze znakiem); jako generator `e^{iKt}` daje poprawną ewolucję unitarną (CTQW). Legacy: 73% ujemnych wag zamyka interpretację stochastyczną; zostaje wyłącznie hamiltonian ze znakowanymi sprzężeniami.
- **Czy któreś przestaje być funkcją Greena?** Oba nigdy nią nie były (Część I). Ale w słabszym sensie: legacy zachowuje przynajmniej *asymptotykę* Greena 3D w nieskończoności; strict porzuca nawet asymptotykę (d^−1,8, „wymiar efektywny 4,6") i w całości przechodzi po stronie operatora.
- **Która wersja jest naturalniejsza z punktu widzenia teorii operatorów?** **Strict.** Samosprzężone, lokalne, elementowo nieujemne, o symbolu eliptycznym pierwszego rzędu z gładkim wycięciem — to podręcznikowy „dobry" operator na grafie: generator dyfuzji/drgań, legalny laplasjanowański obiekt teorii grafów spektralnych i GSP (graph signal processing). Legacy jako operator jest dystrybucyjny w naturze (powłoka spektralna) — obiekt szlachetny, ale inny: to raczej **wzmacniacz rezonansowy** niż operator propagacji.

**Konsekwencja dla pytania badawczego:** przejście legacy → strict w repo (marzec 2026) było — niezależnie od intencji autora — przejściem od *obiektu on-shell typu „widmo/rezonans"* do *obiektu off-shell typu „operator/medium"*. To właśnie strict jest właściwym punktem startowym dla adaptacyjnej interpretacji: operatory się adaptuje, projektory spektralne się wybiera.

---

## CZĘŚĆ III — CZY JĄDRO MOŻE SIĘ UCZYĆ?

### 3.1. Czy konstrukcja K(x,y) → K_t(x,y) jest matematycznie spójna? — TAK, z warunkami

Formalnie: sprzężony układ

```
ψ̇ = i K_t ψ          (szybka dynamika stanu; unitarna, gdy K_t = K_t†)
K̇_t = F(ψ, K_t)      (wolna dynamika operatora; „plastyczność")
```

jest dobrze postawiony jako układ ODE na ℋ × Sym(ℋ), o ile F jest lokalnie lipschitzowska. Struktura operatora zachowana, jeśli F(ψ,·) trzyma przestrzeń symetrycznych (hermitowskich) operatorów — co reguły typu Hebba (F = η(ψψᵀ − γK)) gwarantują. Utrzymywane są: symetria, zerowa diagonala (przez rzut), lokalność (przez maskę grafu lub szybki zanik). **Nie jest** automatycznie zachowana: cyrkulantowość/translacyjność (Hebb ją psuje poza średnią stacjonarną), dodatnia określoność, skala (norma) — stąd obowiązkowy człon zaniku/normalizacji.

Zastrzeżenie fizyczne: czyniąc K zmienną w czasie, `ψ̇ = iK_tψ` traci dokładną unitarność tylko w rzędzie O(K̇) — przy rozdzieleniu skal czasowych (szybki ψ, wolne K) unitarność zachowana adiabatycznie. To standardowa konstrukcja adiabatyczna.

### 3.2. Mapa na istniejącą literaturę (bez ograniczania do fizyki)

| Konstrukcja FIN-adaptacyjna | Najbliższy dojrzały odpowiednik | Komentarz |
|---|---|---|
| K jako operator całkowy na grafie/chmurze punktów | **Neural Operators / Graph Kernel Networks** (Li i in. 2020); **FNO** | warstwa GNO to dokładnie ψ ↦ σ(∫K(x,y)ψ(y)dy); FIN = **jedna warstwa GNO z ręcznie ustawionym, nieuczonym jądrem** |
| Symbol λ(μ) jako odpowiedź częstościowa | **Graph Signal Processing**: K jest *filtrem grafowym* | strict: filtr dolnoprzepustowy z zerem; legacy: pasmowy |
| K̇ = η(ψψᵀ − γK) | **kernel machines z adaptacją jądra**, reguła Oji; w fizyce: **samokonsystentna funkcja Greena** (polaryzacja ośrodka przeubrania propagator — SCBA) | najgłębsza fizyczna analogia: ośrodek, którego propagator sam się przeubrania pod wpływem wzbudzenia |
| ψ̇ = iKψ z K stałym | **reservoir computing** (zbiornik liniowy) / **CTQW** | FIN jest rezerwuarem bez warstwy odczytu |
| K_t adaptowane przez dynamikę | **dynamic graph Laplacians / Graph Structure Learning** | GSL: wspólne uczenie adjacency i dynamiki — dokładnie to samo zadanie |
| Chmura punktów + macierz jądra | **metody Nyströma**, euklidesowe macierze losowe, diffusion maps | FIN-numereryka to Nyström bez świadomości |
| Ciagła wersja K(x,y) | **continuous neural fields / continuous convolutions** | jądro jako funkcja odległości = neural field |
| Wnioskowanie K z trajektorii ψ | **operator inference / system identification** | K jest identyfikowalne z pełnych trajektorii (linear dynamics) — zadanie odwrotne dobrze postawione |
| K(x,y) zależne od stanu ψ | **kernel attention** (Transformery): A(ψ) = softmax(QKᵀ/√d) | patrz Część V |

**Werdykt Części III:** hipoteza adaptacyjności jest spójna i ma **gotową infrastrukturę matematyczną** w co najmniej czterech literaturach (GNO/FNO, GSP, GSL/dynamiczne laplasjanu, samokonsystentne funkcje Greena). Nikt nie musi jej wymyślać — trzeba ją tylko poprawnie zainstancjonować. Repozytorium ma zresztą własny zalążek: pakiet **P2772** testował jawne prawo uaktualnienia `θ_{t+1} = θ_t − lr·∇L_geo(θ)` dla parametrów jądra i uczciwie certyfikował **niestacjonarność** obu krotek parametrów (min ‖grad‖ = 0,087) — ale nie wyciągnął wniosku, że gradient na parametrach to nie jest reguła lokalna typu Hebb i że brak stacjonarności wobec tej jednej funkcji kosztu nic nie mówi o innych.

---

## CZĘŚĆ IV — HEBBIAN LEARNING: ANALIZA WARUNKÓW MATEMATYCZNYCH

Pytanie: czy istnieje naturalna reguła `dK_ij/dt = η ψ_i ψ_j` zachowująca strukturę operatora, stabilna, prowadząca do samoorganizacji, pamięci, emergencji modów, kompresji spektralnej i uczenia reprezentacji? Odpowiadam po punktach, z weryfikacją numeryczną na zamrożonych jądrach.

### 4.1. Struktura i stabilność (zmierzone)

- **Czysty Hebb** (bez zaniku): **niestabilny strukturalnie** — ‖K_t‖ rośnie liniowo bez ograniczeń (numerycznie: 2,5 → 22 w 200 krokach przy η = 0,01). Brak punktu stałego. Odpada.
- **Hebb + zanik Oji** (`K̇ = η(ψψᵀ − γK)`): stabilny, ograniczony; zachowuje symetrię i zerową diagonalę. Ale: dla stacjonarnego pobudzenia zbiega do **projektora na modę dominującą** — zmierzona efektywna ranga po zbiegnięciu: **1,01** (korelacja z ψψᵀ = 1,000). Kompresja spektralna: TAK — do rangi 1. To jest twierdzenie Oji (1982): reguła Hebba z zanikiem realizuje **PCA** — zbiega do głównej składowej kowariancji pobudzenia.
- **Kluczowy fakt o samoodniesieniu:** gdy pobudzenie pochodzi z własnej dynamiki jądra (`ψ(t) = e^{iKt}ψ₀`), średnia czasowa kowariancja to `C = Σ_m |c_m|² P_m` (projektory spektralne K ważone okupacją, fazy się depfazują). Reguła Hebba uczy się **C, nie K**: zmierzone corr(W_nauczone, C) = **0,999**, ale corr(W_nauczone, K_strict) = **0,146**. **Jądro NIE jest punktem stałym uczenia hebbowskiego napędzanego własną dynamiką.** Punktami stałymi są ważone mieszanki jego projektorów spektralnych — to odpowiada na pytanie o „pamięć": system pamięta *statystykę własnych wzbudzeń*, nie *siebie*.
- **Cyrkulantowość:** psuta przez Hebba poza średnią; zachowana w oczekiwaniu przy pobudzeniu stacjonarnym translacyjnie.

### 4.2. Samoorganizacja, pamięć, mody, reprezentacje

- **Samoorganizacja:** TAK — łamanie symetrii przez warunek początkowy/okupację modów (W kolapsuje na kierunek dominujący).
- **Pamięć:** TAK — W jest wykładniczą średnią kroczącą korelacji przeszłych; z energią `E = −½ψᵀKψ` to dokładnie **maszyna Hopfielda**: wzorce zapisane jako mody własne; pojemność ~0,14N (klasycznie).
- **Emergencja NOWYCH modów własnych:** **NIE dla reguł liniowych** (Oja/Sanger potrafią wzmacniać jedynie kierunki obecne w kowariancji pobudzenia — span(C)); TAK dla reguł nieliniowych (BCM z ruchomym progiem → selektywność konkurencyjna; STDP ze strukturą opóźnień → sekwencje/kierunkowość) lub z progowaniem. To ostry warunek: kto chce „nowych modów", musi wyjść poza Hebba liniowego.
- **Kompresja spektralna:** TAK (udowodnione numerycznie do rangi 1; z regułą Sanger/GHA — do kontrolowanej rangi r).
- **Uczenie reprezentacji:** TAK, w sensie PCA/spektralnym (Oja = online PCA); to jest uczenie reprezentacji liniowej.

### 4.3. Krytyczna ocena „hebbowskości" samego repo

Repo zawiera `nadsoliton_neural_analysis.py`, który „dowodzi" emergencji K_legacy z reguły Hebba. Niezależna replikacja jego protokołu (pobudzenie `cos(ωi + t)` + szum, 30 tys. kroków): nauczone W ma corr = **0,9998 z czystą macierzą Toeplitza cos(ω(i−j))**, ale tylko **0,39 z faktycznym K_legacy**. Reguła nauczyła się *pobudzenia*, nie *jądra*; obwiedni tłumiącej (1+βd)⁻¹ ani fazy φ = π/6 **nie da się** wyciągnąć z pobudzenia rezonansowego bez dodatkowej struktury (zmienna β = 1: corr spada do 0,41). Wniosek skryptu jest zaszyty na sztywno w kodzie (`"SIMULATION VERIFIED"`). To nie obala idei — obala konkretny „dowód".

### 4.4. Ranking reguł dla adaptacyjnego K

| Reguła | Stabilność | Co daje | Co psuje | Ocena dla FIN |
|---|---|---|---|---|
| Hebb czysty | ✗ | nic stabilnego | normę | odpada |
| **Oja / Sanger (GHA)** | ✓ | PCA, kompresję do rangi r | nic | dobry fundament; brak nowych modów |
| **BCM** (ruchomy próg θ_M ∝ ⟨ψ²⟩²) | ✓ | selektywność, konkurencję modów, nowe kierunki | symetrię przy braku symetryzacji | najlepszy kandydat „samoregulującego się jądra" |
| STDP | ✓ | sekwencje, kierunkowość, przyczynowość | **symetrię K → niehermitowskość → nieunitarność ψ̇** | wymaga symetryzacji lub rozszczepienia K = K⁺ + K⁻ |
| **Predictive coding / free-energy**: K̇ = −∂F/∂K | ✓ (Lyapunov z definicji) | gwarantowaną stabilność, wariacyjny sens | nic istotnego | **najbardziej zasadnicza**: para (ψ, K) jako wspólny spadek gradientowy; macierz precyzji Λ = (odwrotność kowariancji) jest adaptowaną regułą hebbowską — a **macierz precyzji to właśnie „odwrotna funkcja Greena"** |
| Energy-based (Hopfield nowoczesny, Krotov–Hopfield) | ✓ | pamięć asocjacyjną o dużej pojemności, most do attention | liniowość dynamiki | bardzo mocny kandydat (Część V/VI) |

---

## CZĘŚĆ V — CZY TO JEST UKRYTA SIEĆ NEURONOWA?

Eksperyment myślowy: widzę wyłącznie `ψ̇ = iK_t ψ`, `K̇ = η(ψψᵀ − γK)`, K cyrkulant oscylacyjny na chmurze punktów. Jako artykuł z ML czytam to tak:

| Kandydat | Podobieństwa | Różnice | Werdykt |
|---|---|---|---|
| **Attention / Transformer** | macierz wag na węzłach; mieszanie informacji między pozycjami | attention jest **zależny od danych**: A(ψ) = softmax(Q(ψ)K(ψ)ᵀ/√d) — tzn. Transformer JEST adaptacyjnym jądrem; FIN ma jądro zamrożone, bez softmax, bez content-addressing | **FIN = „zamrożony attention"**; hipoteza adaptacyjna = dodanie zależności K(ψ) — czyli dokładnie powrót do attention |
| **GNN / message passing** | ψ ← Kψ to jedna runda linearny message passing; K_strict elementowo ≥ 0 ⇒ legalne wagi | brak nieliniowości σ, brak wag uczonych per-kanał | jednowarstwowy liniowy GNN |
| **Spectral GNN / ChebNet** | symbol λ(μ) = odpowiedź częstościowa filtru spektralnego | w GNN współczynniki filtru są uczone; tu ustawione | FIN = spectral GNN z ręcznym filtrem |
| **Hopfield / pamięć asocjacyjna** | energia −½ψᵀKψ; mody własne = wzorce | Hopfield ma nieliniową dynamikę retrievalu; FIN liniową | K jest macierzą Hopfielda w reżimie liniowym |
| **Neural fields / continuous kernels** | K(d) funkcją odległości — ciągłe jądro | brak parametryzacji neuronowej | tak, continuous convolution |
| **Continuous attractor** | translacyjna symetria cyklu + lokalność strict ⇒ bump attractor | brak nieliniowości stabilizującej bump | szkielet attractora bez stabilizacji |
| **Graph diffusion / random walk** | strict po normalizacji wierszowej ≈ operator przejścia | znaki (strict ich nie ma — OK; legacy ma 73% ujemnych — wyklucza) | strict: tak; legacy: nie |
| **Kernel PCA / diffusion maps** | wektory własne K = embedding | kernel PCA wymaga PSD; K nie jest PSD | dopiero po PSD-ifikacji (K², \|K\|, e^{tK}) |
| **Gaussian processes** | macierz jądra na chmurze punktów | **nie PSD ⇒ nie jest kowariancją GP** | nie; to jądro nieokreślone (przestrzenie Kreina — niszowa, ale istniejąca literatura) |
| **Nyström** | losowe punkty + jądro = aproksymacja operatora całkowego | — | dokładnie to |
| **Graph transformer** | pozycyjne kodowanie spektralne (wektory własne laplasjanu) | brak uczenia | analogia strukturalna |
| **Reservoir computing** | losowy graf + liniowa dynamika | brak readoutu | rezerwuar liniowy |

**Najtrafniejsza jedna etykieta:** **jednowarstwowy Graph Kernel Network (neural operator) z zamrożonym jądrem i liniową aktywacją** — albo równoważnie **CTQW / rezerwuar liniowy**. „Ukrytą siecią neuronową" FIN jest w tym samym sensie, w jakim każdy operator całkowy jest siecią: architektura jest, uczenia (poza zalążkami z Cz. III–IV) nie ma.

---

## CZĘŚĆ VI — CZY ISTNIEJE NOWA KLASA MODELI?

Szukana klasa: **funkcja Greena + adaptacyjne jądro + graf + lokalna reguła uczenia**.

**Komponenty istnieją i są dojrzałe:** (1) Neural Operators/GNO/FNO — uczą właśnie „funkcji Greena" z danych jako operatorów całkowych na grafach; (2) Graph Structure Learning — adaptacja adjacency; (3) nowoczesne pamięci asocjacyjne (Krotov–Hopfield) — energia + lokalna reguła; (4) predictive coding z uczeniem macierzy precyzji — adaptacyjna „odwrotna funkcja Greena" z lokalną regułą; (5) samokonsystentne funkcje Greena w fizyce (SCBA).

**Konkretna kombinacja FIN — samosprzężone oscylacyjne jądro jako hamiltonian, unitarna dynamika stanu, wolna plastyczność hebbowska z homeostazą — jako nazwana klasa:** nie znalazłem jej jako ustanowionej rodziny. Najbliższe nazwane obiekty:
- **„Hamiltonian neural network with Hebbian plastic coupling"** — opis, nie nazwa;
- **Vanchurin (2020) „world as neural network"** — to samo małżeństwo dynamiki hamiltonowskiej i uczenia, ale przez funkcjonał wariacyjny, nie przez regułę lokalną na jądrze;
- **predictive coding networks z uczoną precyzją** — najbliżej strukturalnie (precyzja ↔ odwrotny propagator, reguły lokalne), ale bez unitarności;
- **equilibrium propagation** (Scellier–Bengio 2017) — energetyczna sieć dynamiczna z quasi-lokalną regułą uczenia.

**Czym różniłaby się klasa FIN-adaptive od znanych:** (i) **oscylacyjne, nieokreślone jądro** jako stan (nikt w ML nie uczy celowo jąder nieokreślonych — przestrzenie Kreina zamiast RKHS); (ii) **unitarna szybka dynamika** zamiast gradientowej/dyfuzyjnej; (iii) **samoodniesienie**: jądro uczy się statystyki własnych wzbudzeń (bootstrap operatorowy `K = F(C(K))`) — to jako zadanie punktu stałego jest dobrze postawione i, o ile wiadomo, niezbadane systematycznie; (iv) hamiltonowska struktura daje prawa zachowania, których modele ML zwykle nie mają. To wystarczy na jedną publikowalną pracę typu „Hebbian dynamics on oscillatory graph operators: fixed points, spectral compression and self-referential bootstrap" — ale nie na „nowy paradygmat".

---

## CZĘŚĆ VII — KONKRETNY WERDYKT (wyłącznie jako matematyk)

**1. Czy hipoteza „adaptacyjnego operatora propagacji informacji" ma sens matematyczny?**
**Tak.** Sprzężony układ (unitarna dynamika stanu + wolna plastyczność operatora z zanikiem/normalizacją) jest dobrze postawiony, zachowuje strukturę operatora i ma ostre, policzalne problemy: punkty stałe `K = F(C(K))`, stabilność, widmo emergentne. Część I pokazuje dodatkowo, że właściwym obiektem adaptacji jest **operator (strict), nie funkcja Greena** — a to porządkuje całą konstrukcję: adaptujesz L, propagatorem jest L⁻¹.

**2. Czy wart rozwijania jako osobny kierunek?**
**Warunkowo tak — jako matematyka/ML, nie jako fizyka.** Oś publikacyjna: spektralna teoria oscylacyjnych euklidesowych operatorów losowych + dynamika hebbowska na nich (punkty stałe, kompresja spektralna, bootstrap). To ma literaturę sąsiednią, ostre problemy i realną publiczność (GSP, neural operators, associative memory). Nie warto rozwijać jako rozszerzenia „teorii pola" — brak wyprowadzenia jądra z czegokolwiek czyni każdą fizyczną interpretację dekoracyjną.

**3. Czy repo zawiera zalążki takiego modelu, nawet jeśli autor nie opisał go tak?**
**Tak — i autor częściowo opisał go tak, tyle że bez kontroli matematycznej:**
- `nadsoliton_neural_analysis.py` — reguła Hebba+Oja na „modach próżni" (wnioski błędne numerycznie, patrz §4.3, ale konstrukcja obecna);
- `QW-540` („HEBBIAN GRAVITY") — jawna reguła `dK_ij/dt = η|ψ_i ψ_j| − decay·K_ij` — to jest dokładnie hipoteza Części III, zapisana w kodzie;
- **P2772** — świadkowie prawa uaktualnienia samouczącego się jądra (gradient na parametrach; certyfikowana niestacjonarność);
- rodzina audytów `hebbian_*` (energy landscape, resonance learning, weight stabilizer) w `fundamental_action_reconstruction/`;
- `p2631_s1581_neural_information_flux_beta_criticality_audit` i pokrewne (strumień informacji neuronowej).
Zalążek jest więc nieobecny „nieświadomie" — jest obecny **jawnie, lecz rozproszony i bez aparatu** (brak twierdzeń o zbieżności, brak punktu stałego, błędne wnioski z korelacji).

**4. Elementy kodu sugerujące nieświadome zbliżanie się do tej klasy modeli:**
(i) używanie K jako **hamiltonianu** `e^{iKt}` (CTQW = linia prosta do GNN); (ii) aktualizacje typu message-passing w QW-593; (iii) audyty ram Fouriera (P3137) — to jest *in spe* analiza filtru grafowego; (iv) gradientowa adaptacja parametrów jądra (P2772 — operator learning w zarodku); (v) ciągłe „spectral tripartition" widma (QW-2118) — nienazwana analiza spektralna kernela; (vi) struktura „oktaw × odległość" = embedding wieloskalowy jak w diffusion maps.

**5. Trzy publikacje jako następny krok:**
1. **Z. Li, N. Kovachki, K. Azizzadenesheli, B. Liu, K. Bhattacharya, A. Stuart, A. Anandkumar — „Neural Operator: Graph Kernel Networks for Partial Differential Equations" (2020).** To jest dokładnie formalizm, w którym żyje K: operator całkowy na grafie jako warstwa, jądro jako uczony obiekt, funkcja Greena jako cel. Pokazuje, co znaczy „uczenie funkcji Greena" zrobione dobrze.
2. **E. Oja — „A simplified neuron model as a principal component analyzer" (J. Math. Biology, 1982).** Kanoniczne twierdzenia o tym, do czego reguły Hebba z zanikiem mogą zbiegać (PCA, kompresja spektralna, ranga 1) — czyli precyzyjne granice analogii hebbowskiej z Części IV.
3. **D. Krotov, J. Hopfield — „Dense Associative Memory for Pattern Recognition" (NeurIPS 2016).** Nowoczesna pamięć asocjacyjna: energia, lokalne reguły, duża pojemność i — kluczowe — udowodniony później związek z attention. To most między „jądrem jako operatorem" a „jądrem jako siecią".
*(Rezerwowi, gdyby szukać głębiej w stronę widma: Goetschy & Skipetrov, „Euclidean random matrices and their applications in physics" (2013); Scellier & Bengio, „Equilibrium Propagation" (2017).)*

---

## ANEKS — REPRODUKOWALNOŚĆ WYNIKÓW NUMERYCZNYCH

Parametry zamrożone z repo: legacy `(α=4ln2, ω=π/4, φ=π/6, β=0,01)` (`RELEASE_6_2`, `nadsoliton_neural_analysis.py`); strict `(ω=0,18575, φ=0,1625, β=1,0, η=1,8)`, cykl Z₁₂, odległość cykliczna, zerowa diagonala (`report_qw2118_ktotal_spectral_tripartition_gate.json`, `QW_2133`). Symbol cyrkulanty: `λ_m = 2Σ_{d=1}^{5}K(d)cos(2πmd/12) + K(6)(−1)^m`. Efektywna ranga: `(Σ|λ|)²/Σλ²`. Dopasowanie `λ_m = c(z−μ_m)e^{−tμ_m}` metodą najmniejszych kwadratów na m = 0…6 (z = 1,179; t = 0,456; c = 1,306; r = 0,993). Testy hebbowskie: 30–60 tys. iteracji, η = 0,02, zanik 0,02, symetryzacja + zerowanie diagonali co krok (protokół identyczny jak w `nadsoliton_neural_analysis.py`).
