# FIN — Release 10.24

# Programy badawcze 255–266

## Miary operatorowe pamięci, identyfikowalność, tomografia prądów i falsyfikacja analogii adaptacyjnych

**Autor:** Krzysztof Żuchowski  
**Afiliacja:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Data:** 28 lipca 2026  
**Język:** polski  
**Licencja:** CC BY 4.0

---

## Konwencja pewności

- **[Proven]** — wynik wynika z jawnego dowodu skończenie wymiarowego albo
  ścisłej tożsamości; obliczenie numeryczne jest tylko certyfikatem.
- **[Strong evidence]** — deterministyczny, odtwarzalny audyt syntetyczny z
  zamrożonym modelem i progiem.
- **[Moderate evidence]** — stabilny wynik w zadanej klasie, lecz zależny od
  dodatkowej definicji lub równoważności.
- **[Speculative]** — skonstruowany program lub obiekt bez wystarczającego
  dowodu.
- **[Refuted]** — jawny kontrprzykład albo no-go w zadanej klasie.

Status matematyczny nie jest statusem fizycznym. W szczególności `[Proven]`
nie oznacza, że operator został zrealizowany przez przyrodę.

# 1. Executive summary

Wykonano wszystkie programy 255–266. Najważniejsze wyniki są następujące.

1. Dla każdego z \(2^{12}-2=4094\) nietrywialnych kontekstów \(E\subset Z_{12}\)
   samoenergia

   \[
   \Sigma_E(z)=A_{EH}(zI+A_{HH})^{-1}A_{HE}
   \]

   jest transformatą dodatniej skończonej miary operatorowej i jest całkowicie
   monotoniczna. **[Proven]**

2. Deklarowany sześciowymiarowy sektor ukryty w podziale even/odd ma tylko
   pięć wymiarów widzialnych przez \(\Sigma_E\). Jeden tryb jest dokładnie
   niewidoczny. Dodanie dowolnych dalszych rozsprzężonych trybów nie zmienia
   samoenergii. **[Proven]**

3. Konteksty i redukcje Schura tworzą kontrawariantną kategorię. Złożenie
   redukcji jest łączne; maksymalny residual w 20 000 łańcuchach wynosi
   \(6.66\times10^{-16}\). Nie udowodniono aksjomatów snopa. **[Proven]**

4. Chiralna podatność pamięci ma zamknięty wzór:

   \[
   \Xi_E
   =B'RB^\ast+BRB'^\ast-BRC'RB^\ast,
   \qquad R=(zI+C)^{-1}.
   \]

   Jest nieparzysta pod inwersją i spełnia jawne ograniczenie normy. Pozostaje
   odbiornikiem zadanego skrętu, nie źródłem orientacji. **[Proven]**

5. Naiwny dynamiczny RG oparty wyłącznie na ścisłym zmniejszaniu kontekstu
   \(12\to6\to3\to1\) nie ma nietrywialnego fixed pointu ani cyklu w kategorii,
   której obiekt zachowuje wymiar. Potrzebna jest nowa, jawna
   size-restoring equivalence. **[Proven]**

6. Zbudowano wieloczasowy bilans rozróżnialności:

   \[
   D_{\rm in}
   =
   \mathcal L_{\rm env}
   +D_{\rm record}
   +D_{\rm conditional}.
   \]

   Bilans teleskopuje z residualem zero, ale nie jest jeszcze termodynamiką.
   **[Proven]**

7. Zamrożona przed testem mapa mostu usuwająca tylko amplitudę
   \(\alpha_{\rm geo}\) nie przenosi legacy do dodatniej klasy generatorów
   strict. Legacy ma \(\lambda_{\min}=-7.7101\), a jego blok ukryty generuje
   biegun na dodatniej osi \(z=5.8248\). Jest to no-go tylko dla
   amplitude-only completion. **[Proven]**

8. Syntetyczny protokół fingerprint-first/calibration-second przechodzi
   300/300 prób przy 50 000 zliczeń na przygotowanie, lecz nie zastępuje
   zewnętrznego pakietu P241 i nie uruchamia P242. **[Strong evidence]**

9. Sam vertex POVM nie identyfikuje prądów. Skonstruowano dwa dodatnie stany
   \(\rho\) i \(\rho^\ast\) o identycznych populacjach i przeciwnych prądach.
   Potrzebne są obserwable interferometryczne. **[Proven]**

10. Atlas false positives pokazuje, że dodatniość Stieltjesa i składanie Schura
    są ogólne dla dodatnich Laplasjanów, a nie swoiste dla FIN. Najbardziej
    swoisty jest fingerprint widmowy. **[Proven]** dla zamrożonych zespołów.

11. Bez niezależnego zegara \((cA,t)\) i \((A,ct)\) są dokładnie
    nierozróżnialne. Pamięć aparatu wymaga danych wieloczasowych, a zmiana
    kształtu generatora jest widoczna w fingerprint. **[Proven]**

12. Benchmark biologiczno-cybernetyczny daje wynik negatywny dla przewagi FIN:
    pojemność pamięci wynosi \(4.0601\) i jest niższa niż we wszystkich 120
    kontrolach, natomiast powrót po zaburzeniu jest najwolniejszy.
    FIN jest odróżnialnym, silnie wygładzającym rezerwuarem, ale nie wykazano
    przewagi obliczeniowej ani biologicznej. **[Strong evidence]**

![Wyniki Stieltjesa i podatności chiralnej](FIN_Programs_255_266_Figures/p255_p258_stieltjes_chiral.png)

# 2. Zamrożony rdzeń i zakres

Pracujemy z rzeczywistym symetrycznym operatorem strict:

\[
W_{xy}=K_{\rm strict}(d(x,y)),\qquad
A=\operatorname{diag}(W\mathbf 1)-W.
\]

Dla:

\[
K_{\rm strict}(d)
=
\frac{\cos(0.18575\,d+0.16250)}{1+d^{1.8}},
\qquad d=1,\ldots,6,
\]

otrzymujemy:

\[
A=A^\ast,\qquad A\mathbf 1=0,\qquad A\succeq0.
\]

Podział \(V=E\sqcup H\) oznacza jedynie operacyjny podział dostępnych i ukrytych
stopni swobody nadsolitona. Nie wprowadza warstwy informacyjnej „pod”
nadsolitonem.

W obliczeniach:

- użyto deterministycznego ziarna `20260729`;
- nie użyto danych laboratoryjnych;
- nie dostrajano progów po wynikach;
- w P261 zamrożono tylko jeden atom mostu;
- zachowano oddzielenie legacy i strict;
- nie użyto Firecrawl.

# 3. Macierz wyników

| Program | Wynik | Status | Najważniejsza granica |
|---|---|---|---|
| P255 | dodatnia miara operatorowa dla 4094 kontekstów | [Proven] | nie dowodzi fizycznego środowiska |
| P256 | minimalny wymiar ukryty \(=5\), a nie 6 | [Proven] | mikroskopowa realizacja nieunikalna |
| P257 | kontrawariantna kategoria Schura | [Proven] | brak twierdzenia o snopie |
| P258 | zamknięty wzór, kowariancja i bound dla \(\Xi_E\) | [Proven] | receiver, nie selector |
| P259 | no-cycle dla \(12\to6\to3\to1\) | [Proven] | brak size-restoring equivalence |
| P260 | wieloczasowy ledger informacji | [Proven] | nie jest energią ani ciepłem |
| P261 | obstruction amplitude-only | [Proven] | nie obala bogatszego mostu |
| P262 | zamrożony audyt fingerprint + skala | [Strong evidence] | wyłącznie dane syntetyczne |
| P263 | vertex-POVM no-go; obserwable prądu | [Proven] | implementacja aparatu zewnętrzna |
| P264 | atlas false positives | [Proven] | zależny od zamrożonych zespołów |
| P265 | iloraz identyfikowalności mechanizmu | [Proven] | brak unikalnego prawa uczenia |
| P266 | matched reservoir benchmark | [Strong evidence] | brak przewagi biologicznej |

# 4. P255 — twierdzenie Stieltjesa dla wszystkich kontekstów

## 4.1. Twierdzenie

Niech \(A\succeq0\) będzie Laplasjanem grafu spójnego, a
\(V=E\sqcup H\), gdzie \(E\neq\varnothing\) oraz \(H\neq\varnothing\).
Oznaczmy:

\[
B=A_{EH},\qquad C=A_{HH}.
\]

Wtedy \(C\succ0\). Dla \(z>0\):

\[
\Sigma_E(z)=B(zI+C)^{-1}B^\ast.
\]

Jeżeli:

\[
C=\sum_{\mu}\mu P_\mu,
\]

to:

\[
\Sigma_E(z)
=
\sum_\mu\frac{\Gamma_\mu}{z+\mu},
\qquad
\Gamma_\mu=BP_\mu B^\ast\succeq0.
\]

Jest to transformata Stieltjesa dodatniej, atomowej miary operatorowej:

\[
\mathsf M_E=\sum_\mu\Gamma_\mu\,\delta_\mu.
\]

Ponadto:

\[
(-1)^n\Sigma_E^{(n)}(z)
=
n!B(zI+C)^{-n-1}B^\ast\succeq0.
\]

**Dowód.** Dodatnia określoność \(C\) wynika z twierdzenia o głównych
minorach zredukowanego Laplasjanu grafu spójnego. Rozkład spektralny daje
reprezentację miarową. Różniczkowanie resolwenty daje wzór na wszystkie
pochodne. Każdy składnik ma postać \(XX^\ast\). \(\square\)

## 4.2. Granice

\[
\lim_{z\to0^+}\Sigma_E(z)=BC^{-1}B^\ast,
\]

\[
\lim_{z\to\infty}z\Sigma_E(z)=BB^\ast.
\]

W audycie even/odd residualy zachowują oczekiwany rząd:

| Granica | Parametr | Residual operatorowy |
|---|---:|---:|
| \(z\Sigma\to BB^\ast\) | \(z=10\) | \(1.44\times10^{-1}\) |
|  | \(z=10^2\) | \(1.59\times10^{-2}\) |
|  | \(z=10^3\) | \(1.60\times10^{-3}\) |
|  | \(z=10^4\) | \(1.61\times10^{-4}\) |
| \(\Sigma\to BC^{-1}B^\ast\) | \(z=10^{-1}\) | \(9.21\times10^{-2}\) |
|  | \(z=10^{-2}\) | \(9.92\times10^{-3}\) |
|  | \(z=10^{-3}\) | \(9.99\times10^{-4}\) |
|  | \(z=10^{-4}\) | \(1.00\times10^{-4}\) |

Audyt wszystkich kontekstów:

| Wielkość | Wynik |
|---|---:|
| liczba kontekstów | 4094 |
| najmniejsza wartość własna bloku \(C\) | 0.1252698327 |
| minimalna wartość własna signed derivatives \(n=0,\ldots,4\) | \(-3.47\times10^{-15}\) |
| maksymalny residual reprezentacji miarowej | \(3.73\times10^{-15}\) |

Ujemna liczba \(10^{-15}\) jest błędem arytmetyki zmiennoprzecinkowej, nie
naruszeniem dodatniości.

# 5. P256 — minimalna realizacja samoenergii

## 5.1. Nowy obiekt: iloraz realizacji ukrytej

Definiujemy:

\[
(C,B)\sim(\widetilde C,\widetilde B)
\quad\Longleftrightarrow\quad
B(zI+C)^{-1}B^\ast
=
\widetilde B(zI+\widetilde C)^{-1}\widetilde B^\ast
\]

dla każdego \(z\) poza biegunami.

Samoenergia identyfikuje:

- widzialne bieguny \(\mu\);
- dodatnie macierze residuów \(\Gamma_\mu\);
- minimalny wymiar

  \[
  n_{\min}=\sum_\mu\operatorname{rank}\Gamma_\mu.
  \]

Nie identyfikuje:

- ortogonalnej bazy sektora ukrytego;
- trybów z \(\Gamma_\mu=0\);
- dowolnych dodatkowych rozsprzężonych trybów.

## 5.2. Wynik FIN even/odd

| Wielkość | Wynik |
|---|---:|
| deklarowany wymiar \(H\) | 6 |
| ranga controllability/observability | 5 |
| minimalny wymiar Stieltjesa | 5 |
| grupy biegunów \(C\) | 4 |
| bieguny widzialne przez \(\Sigma\) | 3 |
| residual rekonstrukcji pole–residue | \(9.55\times10^{-16}\) |
| residual po dodaniu dwóch niewidzialnych trybów | \(4.16\times10^{-17}\) |

Jeden tryb o biegunie około \(1.961406862\) ma residue rank zero. Jest
matematycznie obecny w zadanym bloku \(H\), ale niewidzialny dla wybranego
kontekstu \(E\).

**Wniosek.** Pamięć nie rekonstruuje unikalnego „mikroskopowego środowiska”.
Rekonstruuje jego minimalną klasę wejście–wyjście. **[Proven]**

# 6. P257 — kategoria kontekstów

Definiujemy kategorię \(\mathbf{Ctx}_A(z)\):

- obiekty: niepuste konteksty \(E\subseteq V\);
- dla \(E\subseteq F\) istnieje morfizm redukcji

  \[
  r_{F\to E}:
  \operatorname{Schur}_{V\setminus F}(zI+A)
  \longmapsto
  \operatorname{Schur}_{V\setminus E}(zI+A).
  \]

Kierunek jest przeciwny do inkluzji, dlatego konstrukcja jest
kontrawariantna.

## Twierdzenie o składaniu

Dla \(E\subseteq F\subseteq G\):

\[
r_{F\to E}\circ r_{G\to F}=r_{G\to E}.
\]

Dowód można zapisać na dwa równoważne sposoby:

1. jako łączność eliminacji blokowej Gaussa;
2. jako łączność iterowanego minimum jednej dodatniej formy kwadratowej:

   \[
   \min_{x_{G\setminus E}}Q
   =
   \min_{x_{F\setminus E}}
   \min_{x_{G\setminus F}}Q.
   \]

Certyfikat:

| Test | Wynik |
|---|---:|
| losowe łańcuchy kontekstów | 20 000 |
| maksymalny residual złożenia | \(6.66\times10^{-16}\) |
| residual identyczności | 0 |

Jest to kategoria redukcji. Nie jest jeszcze snopem: nie podano topologii
Grothendiecka, danych pokrycia ani twierdzenia o jednoznacznym sklejaniu.

# 7. P258 — analityczna chiralna podatność pamięci

Niech:

\[
\Sigma(\theta)=B(\theta)R(\theta)B(\theta)^\ast,
\qquad
R(\theta)=(zI+C(\theta))^{-1}.
\]

Ponieważ:

\[
R'=-RC'R,
\]

otrzymujemy:

\[
\boxed{
\Xi
=
B'RB^\ast
+BRB'^\ast
-BRC'RB^\ast
}.
\]

## 7.1. Ograniczenie normy

Jeżeli \(m=\lambda_{\min}(C)>0\), to:

\[
\|R\|\leq\frac1{z+m},
\]

i dlatego:

\[
\|\Xi\|
\leq
\frac{2\|B'\|\|B\|}{z+m}
+
\frac{\|B\|^2\|C'\|}{(z+m)^2}.
\]

Największy zaobserwowany iloraz lewej strony do boundu wynosi
\(0.5596<1\).

## 7.2. Kowariancja

Jeżeli inwersja wysyła \(\theta\mapsto-\theta\), wtedy:

\[
R_E\Xi_E(z)R_E=-\Xi_E(z).
\]

Maksymalny residual na badanej siatce \(z\) wynosi
\(2.68\times10^{-16}\). Residual wzoru analitycznego względem centralnej
różnicy skończonej jest mniejszy niż \(5.43\times10^{-13}\).

To jest chiralna liniowa odpowiedź pamięci. Znak \(\theta\) nadal jest
dostarczony przez rodzinę testową. `QW-2191` pozostaje otwarte.

# 8. P259 — dynamiczny RG i no-go rozmiaru

Zbadano zagnieżdżony ciąg:

\[
Z_{12}
\supset
\{0,2,4,6,8,10\}
\supset
\{0,4,8\}
\supset
\{0\}.
\]

Dla każdego poziomu zachowano:

- wymiar kontekstu;
- znormalizowane kwantyle widma operatora efektywnego;
- udział samoenergii dla \(z=0.1,0.2,0.5,1\).

Odległości kolejnych deskryptorów wynoszą:

\[
0.5589,\qquad0.7420,\qquad1.4113.
\]

## No-go

Jeżeli obiekt RG zawiera wymiar kontekstu, to ścisła redukcja:

\[
12>6>3>1
\]

nie może mieć nietrywialnego punktu stałego ani cyklu. **[Proven]**

Nie jest to no-go dla renormalizacji w ogóle. Brakującym obiektem jest:

\[
\mathcal E_n:
\operatorname{Ctx}_{N_n}\longrightarrow\operatorname{Ctx}_{N_0},
\]

czyli size-restoring embedding/equivalence, wraz z prawem normalizacji
operatora i zmiennej \(z\). Bez \(\mathcal E_n\) słowo „fixed point” porównuje
obiekty różnych typów.

![Przepływ kontekstów i ledger informacji](FIN_Programs_255_266_Figures/p259_p260_rg_information_ledger.png)

# 9. P260 — wieloczasowy tensor bilansu informacji

Niech \(M_k\) będzie wspólnym kanałem dla dwóch hipotez \(p_{k-1},q_{k-1}\):

\[
p_k=M_kp_{k-1},\qquad q_k=M_kq_{k-1}.
\]

Definiujemy stratę kroku:

\[
\ell_k
=
D(p_{k-1}\Vert q_{k-1})
-D(p_k\Vert q_k)
\geq0.
\]

Dla instrumentu zapisującego wynik \(y\) reguła łańcuchowa daje:

\[
D(p_k\Vert q_k)
=
D(r_k^p\Vert r_k^q)
+
\sum_y r_k^p(y)
D(p_k(\cdot|y)\Vert q_k(\cdot|y)).
\]

Stąd nowy obiekt — **Process Information Ledger**:

\[
\boxed{
D_{\rm input}
=
\ell_{k,\rm env}
+D_{k,\rm record}
+D_{k,\rm conditional}
}.
\]

W pięciu krokach:

| Wielkość | Wynik |
|---|---:|
| całkowita utrata rozróżnialności | 1.0483671856 nat |
| suma strat krokowych | 1.0483671856 nat |
| residual teleskopowy | 0 |
| maksymalny residual chain rule | \(6.94\times10^{-17}\) |
| minimalna strata w 1000 dodatkowych parach | 0.0876620 nat |

Ledger rozdziela to, co traci system, co trafia do rekordu aparatu i co
pozostaje w stanie warunkowym. Nie przelicza natów na dżule. Do tego potrzeba
Hamiltonianu, temperatury, kąpieli i protokołu resetu.

# 10. P261 — dynamiczny completion defect

## 10.1. Zamrożona mapa

Zbadano dokładnie jeden atom:

\[
\mathcal C_{\rm amp}:
K_{\rm legacy}
\longmapsto
\frac{K_{\rm legacy}}{\alpha_{\rm geo}}.
\]

Mapa:

- usuwa wyłącznie amplitudę;
- zachowuje legacy
  \(\omega=\pi/4\), \(\phi=\pi/6\), \(\beta_{\rm tors}=0.01\);
- stosuje ten sam przepis \(A=\operatorname{diag}(W\mathbf1)-W\);
- nie dodaje znaku, przesunięcia, dodatniości ani kompresji \(d^\eta\).

## 10.2. Wynik

| Wielkość | Strict | Legacy po \(\mathcal C_{\rm amp}\) |
|---|---:|---:|
| \(\lambda_{\min}(A)\) | \(-2.64\times10^{-16}\) | \(-7.7101445443\) |
| \(\lambda_{\min}(A_{HH})\) | 1.1710910206 | \(-5.8247982640\) |
| dodatni biegun resolwenty | brak | \(z=5.8247982640\) |

Projektowy defect samoenergii:

| \(z\) | \(\|\widehat\Sigma_{\rm strict}-\widehat\Sigma_{\rm legacy}\|_F\) |
|---:|---:|
| 6.5 | 1.2721 |
| 8 | 1.1547 |
| 12 | 1.0242 |
| 20 | 0.9455 |

Wniosek:

\[
\boxed{
\mathcal C_{\rm amp}
\text{ nie jest mostem do dodatniej klasy Stieltjesa strict}
}
\]

**[Proven]** dla tej mapy.

Wynik nie obala mostu zawierającego jawny strict-side sign/shift, zmianę fazy,
częstotliwości, nieliniową kompresję i źródło selektora. Nie rozpoczyna się
role transfer.

![Obstruction mostu i kalibracja](FIN_Programs_255_266_Figures/p261_p262_bridge_fingerprint.png)

# 11. P262 — Calibrated Fingerprint Experiment

Protokół ma kolejność:

1. zamrozić strict fingerprint;
2. zebrać macierz przejścia;
3. przetestować ilorazy widmowe bez jednostki;
4. dopiero później użyć niezależnego zegara do wyboru reprezentanta orbity
   \((cA,\tau/c)\);
5. nie zmieniać modelu po unblindingu.

Audyt syntetyczny:

| Wielkość | Wynik |
|---|---:|
| replikacje | 300 |
| zliczenia na przygotowanie | 50 000 |
| mediana błędu fingerprint | 0.00620 |
| 95 percentyl błędu fingerprint | 0.01050 |
| mediana względnego błędu skali | 0.00325 |
| 95 percentyl względnego błędu skali | 0.00957 |
| joint pass przy progach 0.03/0.02 | 1.000 |

**[Strong evidence]** jako planowanie metody.

Nie dostarczono niezależnie powierzonych zdarzeń. Walidator P241 nie ma czego
przyjąć, więc P242 pozostaje niewykonany.

# 12. P263 — tomografia prądów i pamięci chiralnej

## 12.1. No-go vertex POVM

Skonstruowano:

\[
\rho_\pm
=
\frac{I}{12}
\pm i\epsilon
\left(
|10\rangle\langle2|
-|2\rangle\langle10|
\right),
\qquad \epsilon=0.04.
\]

Oba stany są dodatnie:

\[
\lambda_{\min}(\rho_\pm)=0.0433333.
\]

Mają identyczne populacje:

\[
\operatorname{diag}\rho_+
=
\operatorname{diag}\rho_-,
\]

ale przeciwne prądy czwartej harmonicznej:

\[
C_4(\rho_+)=-0.0150493,
\qquad
C_4(\rho_-)=+0.0150493.
\]

Zatem vertex POVM jest niewystarczający do tomografii prądu. **[Proven]**

## 12.2. Brakujący instrument

Dla \(d=1,\ldots,5\) definiujemy hermitowskie obserwable:

\[
\mathcal J_d
=
i\,d\sum_xW_{x,x+d}
\left(
|x+d\rangle\langle x|
-|x\rangle\langle x+d|
\right).
\]

Wtedy:

\[
C_d(\rho)=\operatorname{Tr}(\rho\mathcal J_d).
\]

Pięć obserwabli ma rangę Grama równą 5 i condition number \(15.17\).
Mogą zostać zmierzone przez interferometryczne bazy fazowe, ale nie przez sam
odczyt wierzchołków.

\(\Xi_E\) wymaga innego eksperymentu: procesowej tomografii generatora przy
dwóch dostarczonych skrętach \(\pm h\). Residual rekonstrukcji centralną
różnicą wynosi \(4.94\times10^{-11}\).

# 13. P264 — atlas false positives

Zamrożono pięć testów:

1. strict fingerprint z tolerancją 0.02;
2. dodatniość pamięci Stieltjesa;
3. składanie kontekstów;
4. scaffold kowariancji inwersyjnej;
5. dodatni kanał kontrakcji informacji.

Wyniki:

| Zespół | Fingerprint | Stieltjes | Schur | Chiral scaffold | Information | Wszystkie |
|---|---:|---:|---:|---:|---:|---:|
| strict target | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| positive circulant | 0.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.00 |
| positive dense | 0.00 | 1.00 | 1.00 | 0.00 | 1.00 | 0.00 |
| signed circulant | 0.00 | 0.24 | 1.00 | 1.00 | 0.03 | 0.00 |
| strict z perturbacją 5% | 0.52 | 1.00 | 1.00 | 1.00 | 1.00 | 0.52 |

![Atlas wyników fałszywie dodatnich](FIN_Programs_255_266_Figures/p264_false_positive_atlas.png)

Wniosek:

- twierdzenie Stieltjesa jest silne, ale nie swoiste dla FIN;
- składanie Schura jest uniwersalną algebrą;
- odbiciowa kowariancja jest wspólna dla cyrkulantów;
- największą swoistość ma pełny fingerprint strict;
- 5% perturbacja przechodzi tylko w 52% prób przy progu 0.02.

Nie wolno przedstawiać własności ogólnych dodatniego Laplasjanu jako
identyfikatora FIN.

# 14. P265 — identyfikowalność Adaptive Memory Geometry

## 14.1. Dokładny no-go skali

\[
\exp[-t(cA)]
=
\exp[-(ct)A].
\]

Bez niezależnego zegara nie istnieje test odróżniający zmianę skali generatora
od zmiany czasu. Residual wynosi dokładnie zero.

## 14.2. Iloraz mechanizmów

Definiujemy **Mechanism Identifiability Quotient**:

\[
\mathfrak Q_{\rm id}
=
\frac{
\{\text{generator, zegar, aparat, pamięć, prawo uczenia}\}
}{
\text{równość wszystkich dopuszczonych rekordów operacyjnych}
}.
\]

Tabela syntetyczna:

| Scenariusz | Drift fingerprint | Defect semigrupy | Wymagany dodatkowy zasób |
|---|---:|---:|---|
| operator statyczny | 0 | 0 | brak |
| skala generatora albo zegara | \(8.88\times10^{-16}\) | 0 | niezależny zegar |
| niejednorodna zmiana generatora | 0.09239 | 0 | fingerprint |
| pamięć aparatu | 0 | 0.09057 | dane wieloczasowe |

Migawki \(A_0,A_1,\ldots,A_r\) nie identyfikują unikalnego pola wektorowego
\(\dot A=F(A,\mathcal R)\). Nieskończenie wiele funkcji \(F\) interpoluje
skończoną trajektorię. Potrzebne są interwencje, rekord aparatu i hold-out.

# 15. P266 — benchmark biologiczno-cybernetyczny

FIN porównano z:

- 60 rezerwuarami o identycznym widmie i losowej orientacji;
- 60 losowymi rezerwuarami symetrycznymi o tej samej liczbie stanów i
  promieniu spektralnym.

Każdy model miał:

- 12 stanów;
- ten sam ciąg wejściowy;
- liniowy readout;
- identyczną procedurę train/test;
- zadanie odtworzenia opóźnień 1–20.

| Metryka FIN | Wartość | Percentyl wśród 120 kontroli |
|---|---:|---:|
| linear memory capacity 1–20 | 4.0601 | 0.000 |
| efektywny wymiar controllability | 1.0698 | 0.0417 |
| kroki powrotu do \(10^{-3}\) | 135 | 1.000 |
| minimalny readout dla 95% Grama | 1 | 0.0917 |

![Benchmark rezerwuaru FIN](FIN_Programs_255_266_Figures/p266_reservoir_benchmark.png)

FIN jest outlierem, ale nie w kierunku przewagi:

- ma najniższą pojemność pamięci w badanym zestawie;
- ma najwolniejszy powrót dla zadanego zaburzenia;
- jego dostępna dynamika wejście–wyjście jest bliska jednowymiarowej.

**[Refuted]** jako hipoteza o automatycznej przewadze nad standardowym
rezerwuarem 12-stanowym w tym benchmarku.

Nie jest obalona możliwość, że inny input, nieliniowy readout albo jawne prawo
adaptacji zmieni wynik. Taki wariant musi jednak pokonać kontrole na danych
hold-out.

# 16. Nowe obiekty teoretyczne O16–O25

| ID | Obiekt | Definicja | Status | Co rozwiązuje |
|---|---|---|---|---|
| O16 | Positive Operator Memory Measure | \(\mathsf M_E=\sum_\mu BP_\mu B^\ast\delta_\mu\) | [Proven] | klasyfikuje dokładną pamięć |
| O17 | Minimal Hidden Realization Quotient | \((C,B)/{\sim_\Sigma}\) | [Proven] | oddziela widzialne środowisko od niewidzialnych trybów |
| O18 | Schur Context Category | \(E\mapsto\operatorname{Schur}_{V\setminus E}(zI+A)\) | [Proven] | typuje poprawne redukcje |
| O19 | Chiral Memory Kubo Tensor | \(\Xi=B'RB^\ast+BRB'^\ast-BRC'RB^\ast\) | [Proven] | łączy flux response z pamięcią |
| O20 | Size-Restoration Obstruction | brak \(\mathcal E_n\) między różnymi rozmiarami | [Proven] | wyjaśnia, dlaczego naiwny fixed point jest źle typowany |
| O21 | Process Information Ledger | \(\ell_{\rm env}+D_{\rm record}+D_{\rm cond}\) | [Proven] | rozdziela utratę i zapis |
| O22 | Amplitude Completion Obstruction | dodatni biegun i ujemne widmo po \(\mathcal C_{\rm amp}\) | [Proven] | zamyka jeden atom mostu |
| O23 | Current Tomography Instrument | \(\{\mathcal J_d\}_{d=1}^5\) + twist process tomography | [Proven]/[Strong evidence] | domyka brak vertex POVM dla prądów |
| O24 | Mechanism Identifiability Quotient | mechanizmy modulo równość rekordów | [Proven] | oddziela zegar, skalę, pamięć i drift |
| O25 | Matched Reservoir Null Benchmark | FIN kontra kontrole o tym samym \(N,\rho\) i protokole | [Strong evidence] | chroni analogie biologiczne przed nadinterpretacją |

Najważniejszym nowym obiektem jest O17. Pokazuje, że „środowisko” odzyskane z
pamięci jest klasą minimalnych realizacji, a nie unikalnym światem
mikroskopowym.

# 17. Próby falsyfikacji

| Teza | Próba zniszczenia | Wynik |
|---|---|---|
| każda \(\Sigma_E\) jest Stieltjesa | wszystkie 4094 konteksty, pochodne do rzędu 4 | przetrwała |
| ukryty sektor jest unikalny | dodano dwa rozsprzężone tryby | obalona |
| diagram kontekstów jest funktorem | 20 000 zagnieżdżonych redukcji | przetrwała |
| \(\Xi\) jest ad hoc | wyprowadzono jako pochodną resolwenty | obalona krytyka ad hoc |
| \(\Xi\) wybiera orientację | inwersja daje parę \(\pm\Xi\) | teza obalona |
| naiwny decimation RG ma fixed point | wymiar \(12>6>3>1\) | teza źle typowana |
| information loss jest energią | brak skali, kąpieli i Hamiltonianu | teza obalona |
| sama amplituda łączy legacy i strict | ujemne widmo i dodatni biegun | teza obalona |
| vertex POVM mierzy prąd | \(\rho,\rho^\ast\) o tych samych populacjach | teza obalona |
| Stieltjes/Schur identyfikuje FIN | dodatnie losowe Laplasjany przechodzą | teza obalona |
| jeden czas rozdziela zegar i generator | \(\exp[-tcA]=\exp[-ctA]\) | niemożliwe |
| FIN ma naturalną przewagę reservoir | 120 matched controls | teza obalona w zadanym benchmarku |

# 18. Relacja do istniejącej matematyki

## 18.1. Mori–Zwanzig i eliminacja stopni swobody

Pamięć po projekcji nie jest nową ideą. O16–O19 umieszczają FIN w skończonym,
dokładnym odpowiedniku formalizmu projekcyjnego. Różnicą jest jawna,
operatorowo-Stieltjesowska klasyfikacja wszystkich kontekstów konkretnego
generatora strict.

## 18.2. Realization theory

O17 jest standardowym typem wyniku teorii realizacji: dane wejście–wyjście
identyfikują minimalną realizację do równoważności, nie dowolne stany
nieobserwowalne. Nowość dla FIN polega na znalezieniu konkretnego niewidzialnego
trybu i na wpisaniu nieidentyfikowalności do granic ontologicznych.

## 18.3. Process tensors

Proces tensorowy uczy, że wieloczasowy proces jest określony przez odpowiedzi na
sekwencje interwencji, a nie przez jedną mapę stanu. O21 i O24 są klasycznym,
skończonym rdzeniem tej samej dyscypliny. Nie są jeszcze kwantowym process
tensorem FIN.

## 18.4. Reservoir computing

Memory capacity jest ograniczona wymiarem rezerwuaru, ale szczegółowa pojemność
zależy od widma, controllability i readoutu. P266 pokazuje, że regularność FIN
nie gwarantuje przewagi. To kontrolowana analogia, nie identyfikacja z mózgiem.

# 19. Rekomendowane programy 267–280

## P267 — Uniqueness of the operator memory measure

Udowodnić jednoznaczność atomowej miary \(\mathsf M_E\) z funkcji
\(\Sigma_E(z)\) oraz stabilność odzyskiwania biegunów przy szumie.

**Prawdopodobieństwo:** 0.92.  
**Stop rule:** nie interpretować biegunów niewidzialnych.

## P268 — Formalny core kategorii Schura

Zapisać O18 w Lean/Mathlib, Coq albo Isabelle: obiekty, morfizmy, identyczność,
łączność i dodatniość.

**Prawdopodobieństwo:** 0.80.

## P269 — Minimalność instrumentu prądowego

Udowodnić minimalną liczbę ustawień POVM potrzebnych do odzyskania pięciu
\(C_d\), uwzględniając straty, dark counts i nieznane fazy.

**Prawdopodobieństwo:** 0.85.

## P270 — Size-restoring RG equivalence

Skonstruować jedną jawną mapę \(\mathcal E_n\), która porównuje różne rozmiary
bez dopasowywania do strict target, albo udowodnić no-go dla lokalnych,
circulant-preserving embeddings.

**Prawdopodobieństwo konstrukcji:** 0.35.  
**Prawdopodobieństwo no-go:** 0.75.

## P271 — Quantum Process Information Ledger

Zastąpić kanały klasyczne instrumentami CP i udowodnić bilans O21 dla relative
entropy procesów z pamięcią.

**Prawdopodobieństwo:** 0.62.

## P272 — Następny pojedynczy atom completion

Po no-go amplitude-only dopuścić dokładnie jeden jawny strict-side atom:
globalny sign/positive shift **albo** nonlinear damping completion. Mapa i
progi muszą być zamrożone przed porównaniem z strict.

**Prawdopodobieństwo bridge:** 0.25.  
**Prawdopodobieństwo wartościowego obstruction:** 0.80.

## P273 — External calibrated fingerprint

Pozyskać pakiet P241 z niezależnym providerem, registrarem i analystą.
Najpierw fingerprint, potem kalibracja, jedna tura P242.

**Prawdopodobieństwo metodologiczne:** 0.75.  
**Prawdopodobieństwo wykonania bez partnera laboratoryjnego:** niskie.

## P274 — Chiral two-flux process tomography

Zaprojektować aparaturę realizującą \(\pm h\), zmierzyć \(\mathcal J_d\) oraz
zrekonstruować \(\Xi_E(z)\) z niepewnością.

**Prawdopodobieństwo wykonalności do oceny przez fizyka:** 0.60.

## P275 — Analityczne false-positive bounds

Wyprowadzić prawdopodobieństwo przejścia strict fingerprint dla zamrożonych
zespołów random Laplacian zamiast polegać tylko na Monte Carlo.

**Prawdopodobieństwo:** 0.72.

## P276 — Two-clock/two-time identifiability theorem

Określić minimalny zestaw dwóch czasów, niezależnego zegara i kontroli, który
oddziela uniform scaling, shape drift i memory defect.

**Prawdopodobieństwo:** 0.90.

## P277 — Causal identification of the adaptive law

Zamrozić skończoną rodzinę praw \(F(A,\mathcal R)\), projekt interwencji i
hold-out. Nie dopuszczać dowolnego universal approximator.

**Prawdopodobieństwo rozstrzygnięcia:** 0.65.

## P278 — Nonlinear matched reservoir benchmark

Powtórzyć P266 dla NARMA, parity, channel equalization i prediction tasks z
jednakowym budżetem parametrów oraz prerejestrowanym rankingiem.

**Prawdopodobieństwo:** 0.80.

## P279 — Conditional thermodynamic bridge

Po dostarczeniu \((H,T,\text{bath},\text{reset protocol})\) zbadać, kiedy O21
staje się bilansem pracy/produkcji entropii. Wynik musi być jawnie
conditioned-on-CA, nie strict.

**Prawdopodobieństwo matematyczne:** 0.75.  
**Prawdopodobieństwo strict source:** niskie.

## P280 — Two-torsor operational section

Połączyć niezależny clock/scale section z dostarczonym orientation resource w
jednym protokole, bez twierdzenia, że jeden zasób generuje drugi.

**Prawdopodobieństwo konstrukcji warunkowej:** 0.85.  
**Strict selector closure:** nie jest celem programu.

# 20. Ranking następnej rundy

| Ranga | Program | Powód |
|---:|---|---|
| 1 | P267 | domyka matematyczną identyfikowalność pamięci |
| 2 | P269 | zamienia no-go vertex POVM w minimalny instrument |
| 3 | P276 | najkrótszy most od ilorazu skali do poprawnego eksperymentu |
| 4 | P268 | formalizuje centralną kategorię |
| 5 | P275 | mierzy rzeczywistą swoistość FIN |
| 6 | P271 | rozszerza ledger do procesów kwantowych |
| 7 | P277 | rozstrzyga, czy adaptacja jest identyfikowalna |
| 8 | P278 | kontynuuje falsyfikację analogii neuronowej |
| 9 | P270 | wysoka nagroda, lecz brak size-restoring map |
| 10 | P272 | tylko jeden nowy atom bridge |
| 11 | P274 | wymaga oceny i aparatury laboratoryjnej |
| 12 | P273 | wymaga rzeczywiście niezależnych danych |
| 13 | P279 | tylko jako teoria warunkowa |
| 14 | P280 | domyka operacyjny pakiet, nie strict selector |

# 21. Granice zgodne z guardrails

Runda 255–266 nie eksportuje:

- niepremisowego selektora;
- rozładowania `QW-2191`;
- kanonicznej jednostki czasu, długości, masy, energii lub działania;
- pełnej mapy legacy→strict;
- role transfer dla legacy;
- źródła \(\beta,\eta,\omega,\phi\) strict;
- unikalnego środowiska mikroskopowego;
- unikalnego prawa adaptacji;
- fizycznego aparatu lub laboratoryjnego rekordu;
- \(L_{\rm total}\), Modelu Standardowego, GR ani ToE.

P261 zamyka wyłącznie amplitude-only completion. Kolejny program mostowy jest
dopuszczalny tylko dla jednego nowego, jawnie typowanego atomu.

# 22. Werdykt końcowy

Najgłębszy wynik tej rundy nie polega na znalezieniu nowej stałej fizycznej.
Polega na zmianie statusu słowa „pamięć”.

\[
\boxed{
\begin{gathered}
\text{Pamięć FIN jest dodatnią miarą operatorową widzianą przez kontekst,}\\
\text{a nie unikalnym ukrytym światem.}
\end{gathered}
}
\]

Samoenergia określa minimalną klasę realizacji wejście–wyjście. Redukcje tych
klas składają się kategoryjnie, chiralny skręt ma dobrze określony tangent,
a informacja posiada dokładny ledger operacyjny. Jednocześnie skala zegara,
źródło orientacji, mikroskopowa realizacja i prawo adaptacji pozostają
nieidentyfikowalne bez dodatkowych operacji.

To wzmacnia FIN jako skończoną teorię operatorową pamięci i procesu. Nie
przekształca jej jeszcze w teorię fizyczną.

# 23. Wybrane źródła porównawcze

1. S. Belyi, E. Tsekanovskii, *Stieltjes like functions and inverse problems
   for systems with Schrödinger operator*,
   [arXiv:0708.0452](https://arxiv.org/abs/0708.0452).
2. Y. T. Lin, Y. Tian, M. Anghel, D. Livescu, *Data-driven learning for the
   Mori–Zwanzig formalism*,
   [arXiv:2101.05873](https://arxiv.org/abs/2101.05873).
3. F. A. Pollock et al., *Operational Markov condition for quantum processes*,
   [arXiv:1801.09811](https://arxiv.org/abs/1801.09811).
4. H. Jaeger, *Short term memory in echo state networks*, GMD Report 152,
   [Fraunhofer record and DOI](https://publica.fraunhofer.de/entities/publication/9dfaead1-4dc0-4e3c-b89b-596f50f671c1).
5. F. Dörfler, F. Bullo, *Kron Reduction of Graphs with Applications to
   Electrical Networks*,
   [arXiv:1102.2950](https://arxiv.org/abs/1102.2950).

Źródła służą do klasyfikacji istniejącej matematyki. Nie są dowodem fizycznej
równoważności FIN z żadnym z tych formalizmów.

# 24. Reprodukcja

```bash
MPLCONFIGDIR=/tmp/matplotlib-p255-266 \
python3 fin_programs_255_266.py

MPLCONFIGDIR=/tmp/matplotlib-p255-266 \
python3 -m unittest -v test_fin_programs_255_266.py
```

Oczekiwany wynik:

- 4094 konteksty P255;
- 20 000 łańcuchów P257;
- 300 replikacji P262;
- 401 wierszy atlasu P264;
- 121 rezerwuarów P266;
- 15/15 testów.
