# P3172 / S2122 — Legacy\* jako generator operatorów i modeli: ocena potencjału

**Data:** 2026-07-20  
**Status:** analysis / potential certificate (nie closure, nie ToE)  
**Artefakty:**  
- `p3172_s2122_legacy_star_operator_model_generator_potential_audit.py`  
- `generated/p3172_s2122_legacy_star_operator_model_generator_potential_audit.json`  
- `generated/p3172_s2122_legacy_star_operator_model_generator_potential_audit.md`  
- `test_p3172_s2122_legacy_star_operator_model_generator_potential_audit.py`  

**Powiązania:** P3171 (faza i dualne dynamiki), P3075–P3082 (Dirichlet/continuum), P3145 (reverse SM/GR layout), P3168 (no-strict-unit).

---

## Guardrails (wiązanie)

1. **Nie zakładamy prawdziwości teorii fizycznej.** Legacy\* jest programem badawczym / obiektem matematycznym.
2. **Brak mostu do Strict.** Nie dopasowujemy parametrów do \(K_{\mathrm{strict}}\). Nie transferujemy ról legacy.
3. **Brak promocji** `L_{\mathrm{total}}`, SM/GR, selektora, ToE z samego \(K^*\).
4. Przyjmujemy wyłącznie fakty robocze z zapytania:
   - to samo jądro generuje \(U(t)=e^{-itL}\) i \(T(t)=e^{-tL}\);
   - nie ma paradoksu obserwatora jako sprzeczności matematycznej: wybór jest modelem operacyjnym;
   - Legacy\* jest samodzielne.

---

## Definicja obiektu

\[
K^*(d)=\frac{A\cos\bigl(\pi d/4+\pi/6\bigr)}{1+\beta d},\qquad A>0,\ \beta>0.
\]

**Realizacje operatorowe (kanoniczne, skończone):**

| Symbol | Definicja | Komentarz |
|--------|-----------|-----------|
| \(W[K^*]\) | \((W\psi)_i=\sum_j K^*(d(i,j))\psi_j\) | undirected radial circulant na nośniku \(V\) |
| \(W^+\) | \(\max\{K^*,0\}\) off-diagonal | rzut dodatni |
| \(L=\mathrm{Deg}-W^+\) | Laplasjan grafowy | \(L\succeq 0\), \(L\mathbf{1}=0\) |
| \(Q=W^+-\mathrm{Deg}\) | generator Markowa | wiersze sumują się do 0 |
| \(\mathcal{A}=sI-W\) | generator przesunięty | \(s\ge\lambda_{\max}(W)\) dla \(\mathcal{A}\succeq 0\) gdy \(W=W^*\) |

Domyślny nośnik obliczeniowy: **cykl \(Z_{12}\)** (skończony, Fourier-diagonalny). Twierdzenia poniżej są **scoped** do tej geometrii, o ile nie napisano inaczej.

---

# CZĘŚĆ I — Klasy operatorów

Werdykt: **TAK** / **NIE** / **WARUNKOWO**.  
Źródło: macierz skończona w JSON + analiza funkcjonalna.

| Klasa | Werdykt | Uzasadnienie |
|-------|---------|--------------|
| **Hermitowski** | **WARUNKOWO** | Undirected \(W[K^*]\) jest rzeczywisty i symetryczny dla wszystkich \(A,\beta>0\) na \(Z_{12}\). Directed residual \(k((j-i)\bmod n)\) jest **niehermitowski**. Hermityzacja jest wyborem geometrii, nie tw. z samego wzoru \(K^*\). |
| **Gram** | **WARUNKOWO** | \(W\) jest jądrem Gram/kowariancji iff \(W\succeq 0\). Dokładnie zabisectionowany próg \(\beta_*\approx 1.07515077623\) na undirected \(Z_{12}\) (A-niezależny, bo \(\mathrm{spec}(AW)=A\,\mathrm{spec}(W)\)). Dla \(\beta<\beta_*\): indefinite Gram. |
| **Laplasjan grafowy** | **TAK** | \(L=\mathrm{Deg}-W^+\) jest zawsze PSD na skończonych próbkach, \(L\mathbf{1}=0\), luka \(\lambda_2(L)>0\) przy spójnym \(W^+\). |
| **Generator Markowa** | **WARUNKOWO** | \(Q=W^+-\mathrm{Deg}\) ma off-diag \(\ge 0\) i wiersze 0. Surowy signed \(K^*\) **nie** jest generatorem Markowa bez \(W^+\). |
| **Schrödinger** | **WARUNKOWO** | \(H=L\) lub \(H=\mathcal{A}\) daje \(U(t)=e^{-itH}\) unitarne przy \(H=H^*\). Residual unitarności numerycznej \(\sim 10^{-15}\). **Nie** eksportuje Born, pomiaru, \(\hbar\). |
| **Dirac** | **NIE** | Brak Clifford, spinorów, pierwiastka z Laplasjanu. Fundamentalny brak struktury. |
| **Maxwell** | **NIE** | Brak 2-formy \(F\), wiązki \(U(1)\), \(dF=0\), \(d{\star}F=J\). |
| **Yang–Mills** | **NIE** | Brak grupy struktury \(G\), połączenia, krzywizny \(F=dA+A\wedge A\). |
| **Lindblad** | **WARUNKOWO** | Mając \(H\), forma Lindblada jest zapisowalna, ale kanały \(V_k\) **nie** wynikają unikalnie z radialnego \(K^*\). |
| **Fokker–Planck** | **WARUNKOWO** | Łańcuch Markowa z \(Q\) jest dyskretnym analogiem FP. Continuum FP wymaga embeddingu i jednostek (P3082: no continuum functor bez importu). |
| **Perron–Frobenius** | **WARUNKOWO** | \(W^+\) nieujemny ⇒ PF przy nieredukowalności. Signed \(K^*\) nie jest operatorem PF. |
| **Koopman** | **WARUNKOWO** | \(U(t)\) jest unitarną ewolucją typu Koopman na \(\ell^2(V)\). Klasyczny przepływ miarowy na fazie **nie** jest rekonstruowany z \(K^*\). |
| **Transport** | **WARUNKOWO** | Directed residual daje circulant transport-like; undirected nie koduje zorientowanego przepływu. |
| **Dyfuzja** | **TAK** | \(T(t)=e^{-tL}\) jest semigrupą dyfuzyjną na \(Z_{12}\) (promień spektralny \(\le 1\), rachunek spektralny dokładny). |
| **Falowy** | **WARUNKOWO** | Równanie falowe \(\partial_t^2 u=-Lu\) wymaga liftu fazowego \((u,\partial_t u)\); nie wynika z samej semigrupy pierwszego rzędu. |
| **Całkowy** | **TAK** | \(W[K^*]\) jest z definicji operatorem integralnym/macierzowym z jądrem \(K^*\). |
| **Pseudoróżniczkowy** | **WARUNKOWO** | Na \(Z_{12}\) circulant = mnożnik Fouriera (ΨDO rzędu 0 na okręgu). Na \(\mathbb{R}^n\) wymaga continuum. |
| **Spektralny** | **TAK** | Tw. spektralne skończone: \(W,L\) diagonalizowalne; rachunek funkcyjny daje \(U,T\). |

### Podsumowanie I (skończone)

- **TAK:** 4 — Laplasjan, dyfuzja, całkowy, spektralny.  
- **WARUNKOWO:** 11.  
- **NIE:** 3 — Dirac, Maxwell, Yang–Mills.

**Wniosek I:** Legacy\* jest **silnym generatorem skończonej teorii spektralnej / grafowej / Markowa / dyfuzji**. Nie jest generatorem teorii cechowania ani spinu.

---

# CZĘŚĆ II — Dualność dynamik

### Założenie robocze

\[
U(t)=e^{-itL},\qquad T(t)=e^{-tL}.
\]

### Czy to dwie reprezentacje jednej głębszej struktury?

**TAK (matematycznie, scoped).**  
Jeśli \(L=L^*\) na skończonej przestrzeni Hilberta, istnieje wspólna spektralna miara \(E\):

\[
U(t)=\int e^{-it\lambda}\,dE(\lambda),\qquad
T(t)=\int e^{-t\lambda}\,dE(\lambda).
\]

Obie dynamiki to **funkcje Borela** tego samego generatora. Głębsza struktura: **rachunek funkcyjny** \(f\mapsto f(L)\).

### Wspólna algebra?

**TAK:** \(C^*(L)\) (przemienna \(C^*\)-algebra generowana przez \(L\)) lub algebra von Neumanna \(\{L\}''\).  
Na \(Z_{12}\) undirected: diagonalizacja Fourierowska — algebra jest izomorficzna z funkcjami na spektrum.

### Wspólny generator uniwersalny?

**TAK:** \(L\) (lub dowolny afiniczny \(sI-W\) o tych samych projektorach spektralnych).  
Uniwersalność: wszystkie dynamiki postaci \(f(L)\) dla Borelowskiego \(f\).

### Klasyfikacja dynamik z jednego operatora

| Typ | Warunek na \(f\) |
|-----|------------------|
| Grupa unitarna | \(f(\lambda)=e^{-itg(\lambda)}\), \(g\) rzeczywiste na \(\mathrm{spec}(L)\) |
| Semigrupa kontrakcji | \(\mathrm{Re}\,f\le 0\) w sensie generatora |
| Semigrupa Markowa | dodatniość + zachowanie prawdopodobieństwa |
| Dyfuzja gradientowa | \(f(\lambda)=-t\lambda\) na PSD Laplasjanie |
| Oscylacja spektralna | \(f(\lambda)=-it\lambda\) |

**Paradoks obserwatora:**  
Nie jest sprzecznością operatorową. Wybór \(U\) vs \(T\) = wybór \(f\) = **model operacyjny** (instrument / filtr / rejestr).  
To jest **postulat modelowania**, nie twierdzenie o świadomości i nie Born-rule.

**Świadek numeryczny (β=1, t=0.25):** unitarity residual \(\sim 10^{-15}\); \(T\) i \(U\) odtwarzane z rachunku spektralnego z residualem maszynowym.

---

# CZĘŚĆ III — Generator modeli (kategoria / funktor)

### Pipeline

\[
\mathrm{Legacy}^* \;\longrightarrow\; \mathrm{Operator} \;\longrightarrow\; \mathrm{Semigroup/Group} \;\longrightarrow\; \mathrm{Model\ fizyczny}.
\]

### Czy to funktor?

**WARUNKOWO.**

- Po ustaleniu kategorii nośników \(\mathcal{C}_{\mathrm{supp}}\) (np. skończone przestrzenie metryczne / cykle),  
  funktora realizacji \(R:\{\text{jądra radialne}\}\to\mathrm{Op}\) (np. \(W\), \(L\), \(Q\)),  
  oraz rachunku \(F_f:\mathrm{Op}\to\mathrm{Dyn}\),  
  złożenie \(F_f\circ R\) **jest** funktorem na tej ustalonej kategorii.
- **Bez** ustalonego nośnika i realizacji Legacy\* jest **generatorem kandydatów**, nie unikalnym funktorem.

### Kategoria operatorów generowanych przez Legacy\*?

**TAK, definiowalna:**

Obiekty (przykłady):  
`W_signed_circulant`, `W_PSD_Gram`, `L_graph_Laplacian`, `Q_Markov`, `H_Schrodinger`, `T_diffusion`, `U_unitary`.

Morfizmy: sploty radialne, ograniczenia do podgrafów, sprzężenia unitarne zachowujące radialność, mapy Borela \(f(L)\).

### Klasa uniwersalności?

**Częściowa:**  
*radial circulant Fourier multipliers with 8-periodic phase modulation and linear damping*.  
To **nie** jest klasa Wigner/RMT (spektrum całkowicie integracyjne na cyklu).

### Klasyfikacja obrazów funktora?

**Częściowa:** w klasie skończonych cykli + ustalonych realizacji — tak.  
Na całe \(B(\mathcal{H})\) — **nie** (obraz jest cienki).

**Nie jest** funktorem ToE.

---

# CZĘŚĆ IV — Odwracalność

### Problem odwrotny: \( \mathrm{Op}\mapsto K^* \)?

| Własność | Werdykt | Komentarz |
|----------|---------|-----------|
| Injektywność na klasie undirected radial circulant \((A,\beta)\) | **TAK (num.)** | różne \((A,\beta)\) dają różne \(W\); profil \(d=0..6\) odzyskiwany z pierwszego wiersza z błędem \(0\) |
| Surjektywność na wszystkie hermitowskie | **NIE** | tylko radialne circulanty o 8-okresowej modulacji i liniowym tłumieniu |
| Odwracalność globalna | **NIE** | brak bijekcji \(\mathrm{Op}\leftrightarrow K^*\) |
| Częściowa odwracalność | **TAK** | na klasie undirected radial \(Z_n\) |

### Porównania ideowe

| Analogia | Podobieństwo | Różnica |
|----------|--------------|---------|
| **Inverse spectral problems** | z widma circulant można odzyskać mnożniki Fouriera, stąd profil radialny | tu znamy całą macierz, nie tylko spektrum; problem jest łatwiejszy w klasie |
| **Gel'fand duality** | dualność spektrum ↔ algebra przemienna | \(C^*(L)\) dualne do \(\mathrm{spec}(L)\), ale **nie** dualność odzyskująca \(K^*\) spoza klasy |
| **Jacobian Conjecture** (tylko idea) | pytanie o globalną odwracalność wielomianową | tu: odwracalność mapy jądro→operator w nieskończenie wymiarowej „przestrzeni jąder”; analogia **ideowa**, nie twierdzenie |

**Wniosek IV:** mapowanie jest **częściowo odwracalne** w naturalnej klasie radialnej; **nie** jest uniwersalnym izomorfizmem.

---

# CZĘŚĆ V — Niezmienniki

### Dynamika unitarna \(U(t)=e^{-itL}\)

| Wielkość | Zachowanie |
|----------|------------|
| Widmo \(L\) | niezmiennicze |
| Norma \(\|\psi\|_2\) | zachowana |
| Klasy sprzężenia unitarnego \(L\) | zachowane w ewolucji operatorowej |
| Energia \(\langle L\rangle\) | zachowana (gdy \(H=L\)) |
| Faza względna modów | ewoluuje; globalna faza swobodna |
| Dodatniość stanów | **nie** (amplitudy zespolone) |
| Topologia grafu | zakodowana w \(L\); nie zmienia się w czasie |
| Informacja von Neumanna czystego stanu | 0; mieszanych: unitarnie niezmiennicza |

### Dynamika dyfuzyjna \(T(t)=e^{-tL}\)

| Wielkość | Zachowanie |
|----------|------------|
| Widmo \(L\) | niezmiennicze (generator) |
| Masa \(\sum\psi\) (dla Markowa) | zachowana przy \(Q\) |
| Dodatniość | zachowana dla semigrupy Markowa z \(Q\) |
| Energia Dirichleta \(\langle L\psi,\psi\rangle\) | maleje (Lyapunov) |
| Entropia względna do stacjonarnego | często maleje (skończone Markow) |
| Okresowość fazy 8 w \(K^*\) | własność jądra, nie dynamiki w czasie |
| Orientacja / selektor | **nie** generowana z undirected \(K^*\) |

### Niezmienniki rodziny \((A,\beta)\)

| Niezmiennik | Zależność |
|-------------|-----------|
| Okres fazy \(\Delta d=8\) | niezależny od \((A,\beta)\) |
| Algebraiczność \(\cos\theta(d)\in\mathbb{Q}(\sqrt{2},\sqrt{3})\) | niezależna |
| Degeneracje par Fourierowskich undirected | tak, dopóki radial |
| Znak spektrum \(W\) | zmienia się przy \(\beta_*\) |
| \(\beta_*\) | niezależne od \(A\) |

---

# CZĘŚĆ VI — Geometria

| Geometria | Status z Legacy\* | Komentarz |
|-----------|-------------------|-----------|
| **Spektralna** | **częściowa TAK** | spektrum \(W,L\); odległości spektralne; brak continuum spectral triple bez importu |
| **Informacji** \(d_I=-\log\|K\|\) | **WARUNKOWO** | **nie metryka** przy małym \(\beta\) (β=0.01: 216/1320 naruszeń trójkąta); **metryka na \(Z_{12}\)** przy β=1 i β=2 (skończony świadek) |
| **Nieprzemienna** | **słaba** | \(C^*(L)\) przemienna na skończonym cyklu; NCG wymaga bogatszej algebry |
| **Grafowa** | **TAK** | ważony cykl, Laplasjan, cięcia, \(W^+\) |
| **Hiperboliczna** | **intuicja / scaffold** | ogon \(1/(1+\beta d)\) przypomina zasięg/hiperbolę w \(d\); **nie** dowodzi geometrii H^n |
| **Emergentna** | **hipoteza** | continuum / Lorentz **nie** wynikają (P3079–P3083 no-go bez importu) |

**Twierdzenie scoped (dI):**  
Definicja \(d_I=-\log|K^*|\) **nie** jest uniwersalnie metryką na rodzinie Legacy\*; reżim dużego \(\beta\) może być metryczny na skończonym nośniku, reżim małego \(\beta\) — nie.

---

# CZĘŚĆ VII — Fizyka

Dla każdej teorii: (a) z samego \(K^*\), (b) z dualności operatorowej, (c) minimalne nowe aksjomaty.

| Teoria | Z Legacy\* | Z dualności \(U/T\) | Minimalne nowe aksjomaty | Poziom eksportu |
|--------|------------|---------------------|---------------------------|-----------------|
| Mechanika klasyczna | profil radialny / jądro | formalne \(H=L\) po imporcie \((q,p)\) | przestrzeń fazowa, forma symplektyczna, jednostka czasu | scaffold |
| Mechanika kwantowa | \(H=H^*\) na \(\ell^2(Z_{12})\) | \(U(t)\) ewolucja | Born, observabla, pomiar, \(\hbar\) | partial_formal |
| Procesy Markowa | \(Q\) z \(W^+\) | \(T(t)=e^{tQ}\) | interpretacja probabilistyczna | **strong_finite** |
| Termodynamika | brak unit-bearing | Lyapunov/entropia formalna | \(T\), kąpiel, miara równowagi | formal_only |
| Teoria informacji | korelacje / \(d_I\) | kanał Markowa | źródło miary prob. | partial |
| Teoria pola | brak QFT | free Laplacian na grafie | continuum, pola, akcja, RG | **absent** |
| Elektrodynamika | brak | brak | \(U(1)\), metryka, prąd | **absent** |
| Standard Model | brak | brak | grupa, reprezentacje, Yukawa, Higgs, jednostki | **absent** |
| OTW | brak | graf ≠ Lorentz | metryka Lorentza, EH, \(T_{\mu\nu}\) | **absent** |
| QG | brak | brak | QFT+GR+Planck | **absent** |

**Wniosek VII:** najsilniejszy eksport fizyczny to **skończone procesy Markowa / dyfuzja spektralna / formalna QM na grafie**. Reszta wymaga dużego stosu aksjomatów **nie wynikających** z \(K^*\).

---

# CZĘŚĆ VIII — Program badawczy (nowe działy?)

Jeśli Legacy\* traktować jako obiekt fundamentalny **matematycznie** (nie ontologicznie), naturalne gałęzie:

| Gałąź | Czy powstaje? | Nazwa robocza |
|-------|---------------|---------------|
| Nowa klasa operatorów | **TAK** | `LegacyStarRadialCirculantFamily` |
| Nowa klasa grafów | **TAK** | `PhaseModulatedWeightedCycles` |
| Nowa klasa jąder | **TAK** | `LegacyStarKernel` |
| Nowa klasa semigrup | **WARUNKOWO** | dualne semigrupy Borela \(f(L)\) |
| Nowa klasa geometrii | **WARUNKOWO** | signed-kernel / reżim-β geometry |
| Nowe problemy odwrotne | **TAK** | radial kernel recovery + β\* threshold |
| Nowa teoria spektralna | **WARUNKOWO** | przejście indefinite→PSD Gram |

**Program czysto matematyczny (rekomendowany):**

1. Analityczne \(\beta_*(n)\) dla cykli \(Z_n\).  
2. Indefinite Gram i Krein-space dynamics dla \(\beta<\beta_*\).  
3. Klasyfikacja \(f(L)\) zachowujących dodatniość / Markowa.  
4. Inverse: charakterizacja profili o okresie fazy 8 i liniowym tłumieniu.  
5. Granica continuum na \(S^1\) jako ΨDO (osobny funktor).

**Program fizyczny** wymaga **jawnego** importu aksjomatów (jednostki, continuum, gauge/spin) — nie „wynikania z \(K^*\)”.

---

# CZĘŚĆ IX — No-go (obalanie programu jako ToE / pełnej fizyki)

| Twierdzenie / cel | Status | Powód formalny | Brakujący aksjomat | Brakująca struktura | Fundamentalne? |
|-------------------|--------|----------------|--------------------|--------------------|----------------|
| Unikalny Dirac z \(K^*\) | **NO-GO** | brak Clifford | spin structure | spinor bundle | **tak** |
| Maxwell / YM z \(K^*\) | **NO-GO** | brak połączenia | zasada cechowania | \(G\)-bundle | **tak** |
| \(d_I\) metryka dla wszystkich β | **NO-GO** | trójkąt pada przy małym β | — | reżim geometryczny | **tak w reżimie** |
| Unikalne kanały Lindblada | **NO-GO** | nieoznaczoność Kraus | system–kąpiel | \(V_k\) | **tak** |
| Continuum QFT | **NO-GO bez importu** | brak functora continuum | embedding + skala | rozmaitość | fund. jako eksport |
| Jednostki fizyczne z \(K^*\) | **NO-GO** | \(K^*\) bezwymiarowe; orbita \((A,\beta)\) | \(S_+/\Omega_{\mathrm{scale}}\) | triada jednostek | **tak na artefaktach** |
| Selektor orientacji z undirected \(K^*\) | **NO-GO** | inversion-even | źródło odd | section torsora | **tak w klasie** |
| SM / GR / ToE z Legacy\* | **NO-GO** | tylko scaffold | duży stos | pełna teoria | **tak** |
| Globalny inverse Op→\(K^*\) | **NO-GO** | niesurjektywność | restrykcja klasy | geometria nośnika | fund. globalnie |
| „Paradoks obserwatora rozwiązany ontologicznie” | **AXIOM** | dualność ≠ teoria świadomości | mapa instrument→\(f\) | POVM | wybór modelowania |

**Wniosek IX:** jako **program unifikacji fizyki** Legacy\* upada bez masywnych dodatkowych aksjomatów.  
Jako **program teorii operatorów / grafów / semigrup** — stoi.

---

# CZĘŚĆ X — Ocena potencjału (0–10)

Ocena = potencjał **generatora badań**, nie prawdziwość ontologiczna.

| Domena | Score | Uzasadnienie |
|--------|------:|--------------|
| Teoria operatorów | **8** | bogata rodzina circulant / dualny rachunek funkcyjny / reżim PSD |
| Analiza funkcjonalna | **7** | semigrupy, grupy, Gram indefinite, progi |
| Geometria spektralna | **6** | grafowa realna; continuum nie |
| Teoria grafów | **7** | ważone cykle, Laplasjan, PF, reżimy |
| Fizyka matematyczna | **6** | czysta dualność \(U/T\); ograniczona fizyka continuum |
| Fundamenty QM | **5** | unitarność tak; Born/pomiar nie |
| Teoria informacji | **5** | korelacje; \(d_I\) reżimowe; kanały Markowa |
| Program unifikacji | **2** | sam w sobie scaffold; SM/GR/QG absent |

**Średnia (arytmetyczna skryptu):** \(\approx 5.75/10\).

### Werdykt końcowy

| Pytanie | Odpowiedź |
|---------|-----------|
| Czy Legacy\* ma potencjał matematyczny? | **Tak — wysoki w spektrum/grafach/semigrupach** |
| Czy Legacy\* ma potencjał jako ToE? | **Nie — bez dużego stosu zewnętrznych aksjomatów** |
| Czy dualność \(U/T\) jest głęboka? | **Tak jako rachunek funkcyjny jednego generatora** |
| Czy rozwiązuje fundamenty fizyki? | **Nie eksportuje jednostek, gauge, spinu, SM, GR** |
| Czy warto rozwijać jako obiekt? | **Tak — jako klasa operatorów, nie jako domknięta fizyka** |

---

## Separacja poziomów wiedzy (obowiązkowa)

### ✓ Twierdzenia udowodnione (scoped, skończone / standard FA)

1. Undirected \(W[K^*]\) na \(Z_{12}\) jest rzeczywisty, symetryczny.  
2. Przy \(L=L^*\): \(U(t)\) i \(T(t)\) mają wspólną miarę spektralną (rachunek funkcyjny).  
3. \(L=\mathrm{Deg}-W^+\) jest PSD, \(L\mathbf{1}=0\).  
4. \(d_I\) **nie** jest metryką przy β=0.01 (216 naruszeń); **jest** metryką przy β=1,2 na \(Z_{12}\) (skończony świadek).  
5. Profil radialny odzyskiwalny z undirected circulant (partial inverse).  
6. Faza \(\cos(\pi d/4+\pi/6)\) ma okres 8 w \(d\).  
7. Dirac / Maxwell / YM **nie** wynikają z samego radialnego skalarnego \(K^*\).

### ✓ Wyniki z analizy (nie pełne tw. continuum)

1. Macierz klas: 4 TAK / 11 WARUNKOWO / 3 NIE.  
2. Legacy\* = silny generator modeli spektralno-Markowskich.  
3. Nie jest unikalnym funktorem fizycznym ToE.  
4. „Brak paradoksu obserwatora” = fakt o wspólnym generatorze + postulat operacyjny.

### ✓ Hipotezy

1. Klasa 8-okresowych damped radial multipliers ma sensowny continuum ΨDO limit na \(S^1\).  
2. Reżim indefinite Gram jest użyteczny poza kowariancją PSD.

### ✓ Intuicje

1. Instrument wybiera \(f\) w \(f(L)\) jak filtr spektralny.  
2. \(1/(1+\beta d)\) = scaffold zasięgu; cosinus = scaffold interferencji.

### ✓ Spekulacje

1. Jakakolwiek głębsza unifikacja „z Legacy\*” wymaga **nowych** aksjomatów źródłowych spoza \(K^*\).

---

## Next honest step

1. **Czysta matematyka:** rozwijać `LegacyStarRadialCirculantFamily` (analityczne \(\beta_*(n)\), Krein, inverse, klasyfikacja \(f(L)\)).  
2. **Fizyka:** tylko z **jawnym** pakietem aksjomatów (jednostki / continuum / gauge-spin), oznaczonym jako import — bez twierdzenia, że wynika z \(K^*\).  
3. **Nie** most Strict, **nie** role transfer, **nie** ToE z dualności \(U/T\).

---

## Reprodukcja

```bash
cd fundamental_action_reconstruction
python3 p3172_s2122_legacy_star_operator_model_generator_potential_audit.py
python3 -m unittest test_p3172_s2122_legacy_star_operator_model_generator_potential_audit.py
```
