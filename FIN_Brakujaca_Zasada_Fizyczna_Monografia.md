# Od operatora widmowego do fizyki?

## Rekonstrukcja brakującej zasady FIN, twierdzenia niemożliwości i minimalny most operacyjny

**Data:** 18 lipca 2026  
**Zakres:** niezależna matematyka i fizyka matematyczna; bez tez kosmologicznych, filozoficznych i „Theory of Everything”  
**Cel:** znaleźć najmniejsze twierdzenie prowadzące od operatora spektralnego do informacji, dynamiki, wymiarów fizycznych i obserwabli — albo udowodnić, że takie twierdzenie nie może istnieć bez nowych aksjomatów.

---

## Konwencja pewności

Każdy akapit zawierający wniosek merytoryczny dziedziczy najbliższą z poniższych etykiet. Etykieta „udowodnione” oznacza dowód w jawnie podanym modelu i zakresie, a nie prawdziwość bez żadnych hipotez.

- **[Udowodnione]** — wynika z podanego dowodu, standardowego twierdzenia albo dokładnego rachunku skończeniewymiarowego.
- **[Silne dowody]** — wspiera to kilka niezależnych argumentów, lecz nie podano pełnego twierdzenia obejmującego wszystkie możliwe rozszerzenia.
- **[Umiarkowane dowody]** — rozsądna hipoteza badawcza z istotnymi otwartymi warunkami.
- **[Spekulatywne]** — kierunek wart testu, ale bez obecnego dowodu.

Słowo „nowe” oznacza dalej „nowo sformułowane lub zsyntetyzowane w niniejszej analizie”. Nie jest roszczeniem pierwszeństwa w światowej literaturze.

---

# 1. Streszczenie wykonawcze

**[Udowodnione]** Najmniejszym rdzeniem generującym oba jądra, miarę spektralną i rodzinę Greena nie jest samo jądro FIN, lecz **dwuoznaczony samosprzężony układ spektralny**

\[
\Sigma_0=(\mathcal H,L,F),\qquad
F:\{\mathrm{s},\mathrm{l}\}\times\sigma(L)\longrightarrow\mathbb R,
\]

gdzie

\[
K_{\mathrm{strict}}=F(\mathrm{s},L),\qquad
K_{\mathrm{legacy}}=F(\mathrm{l},L),\qquad
G_L(z)=(L-zI)^{-1}.
\]

Miara projekcyjna \(E_L\) jest dana automatycznie przez twierdzenie spektralne. Aby rdzeń generował także prawo adaptacji, trzeba dodać pole wektorowe \(V\) na dopuszczalnej przestrzeni operatorów:

\[
\Sigma=(\mathcal H,L,F;\mathcal M,V),\qquad \dot K=V(K).
\]

Samo \(L\) nie wyznacza ani dwóch oznaczeń \(F(\mathrm{s},\cdot),F(\mathrm{l},\cdot)\), ani \(V\).

**[Udowodnione]** Na cyklu \(C_{12}\) oba rozpatrywane radialne jądra są dokładnie wielomianami stopnia co najwyżej sześć w laplasjanie cyklu. Jest to ścisła wspólna struktura spektralna. Nie dowodzi ona jednak unikalnej „głębszej substancji”, ponieważ **każde** rzeczywiste radialne jądro na \(C_{12}\) ma tę własność. Wspólna algebra jest faktem; jej interpretacja jako jedynej przyczyny nie jest wymuszona.

**[Udowodnione]** Brakujący selektor nie może być funkcją obecnych radialnych danych. Inwersja \(R:x\mapsto-x\) pozostawia \(L\), oba jądra, ich resolwenty, funkcje cieplne, zeta i miary spektralne bez zmian, lecz odwraca orientację. Ekwariantna mapa z obiektu ustalonego przez \(R\) do dwuelementowego torsora orientacji nie istnieje. Bifurkacja może wytworzyć parę \(\pm\), ale nie wybiera jednego znaku.

**[Udowodnione]** Nie istnieje naturalna konstrukcja niezerowej wielkości wymiarowej z danych bezwymiarowych. Grupa zmiany jednostek działa trywialnie na wejściu, a nietrywialnie na długości, czasie, działaniu, energii i masie. Ekwariantna mapa do wielkości o niezerowym ciężarze wymiarowym musi być zerowa. Heat kernel, zeta, spectral action, modular flow i RG mogą tworzyć skale **względne**, gdy otrzymają skalę odniesienia, stan albo warunek brzegowy; nie tworzą absolutnej kalibracji eksperymentalnej.

**[Udowodnione]** Informacja Shannona nie wyznacza energii ani temperatury. Dla tego samego rozkładu \(p_i\) i każdego \(\theta>0\) można zbudować

\[
E_i=-\theta\log p_i+C,
\]

co daje ten sam rozkład Gibbsa. Zasada Landauera tworzy poprawny most do ciepła dopiero po zadaniu fizycznego rezerwuaru, Hamiltonianu i temperatury. Jaynes, KMS, względna entropia i termodynamiczne teorie zasobów mają ten sam warunkowy charakter.

**[Udowodnione]** Nie ma pojedynczego „brakującego twierdzenia” zależnego wyłącznie od obecnego operatora, które równocześnie wytworzyłoby selektor, niezerową skalę oraz obserwable. Istnieją co najmniej dwa niezależne twierdzenia no-go: torsor orientacji nie ma punktu stałego, a torsor kalibracji jednostek nie ma kanonicznego elementu. Niejednoznaczność stanów, instrumentów, dynamiki, statystyki cząstek i limitu continuum dodaje trzeci, operacyjny poziom niedookreślenia.

**[Silne dowody]** Najkrótszym realistycznym rozszerzeniem nie jest nowa funkcja spektralna, lecz **punktowana, skalibrowana realizacja operacyjna**

\[
\mathfrak P=
\bigl(
\mathcal A,\mathcal H,K;\,
\xi,\,
\ell_*,t_*,\mathcal A_*,\,
\alpha_\tau,\,
\{\rho_x\},\,
\{\mathcal I_y\},\,
\mathsf{Emb}
\bigr).
\]

Stan lub warunek \(\xi\) może złamać symetrię; \((\ell_*,t_*,\mathcal A_*)\) kalibruje długość, czas i działanie; \(\alpha_\tau\) wybiera dynamikę; \(\rho_x\) i instrumenty \(\mathcal I_y\) definiują przygotowanie i pomiar; \(\mathsf{Emb}\) dostarcza lokalności albo rodziny refinement. To może być aksjomat prawidłowej teorii fizycznej, ale jego składników nie wyprowadza obecny operator.

**[Udowodnione]** Ostateczny wynik misji jest zatem negatywny wobec idei jednego bezzałożeniowego mostu, a pozytywny wobec dokładnie określonego programu: FIN może stać się modelem fizycznym dopiero po podaniu fizycznej realizacji operacyjnej i wygenerowaniu przewidywań bezwymiarowych, których nie użyto do kalibracji.

---

# 2. Co znaczy „stać się fizyką”

**[Udowodnione]** Macierz, widmo, funkcjonał wariacyjny albo równanie ewolucji są obiektami matematycznymi. Aby tworzyły teorię fizyczną, potrzebne są co najmniej:

1. klasa procedur przygotowania;
2. prawo dynamiki z fizyczną parametryzacją;
3. klasa zdarzeń lub instrumentów pomiarowych;
4. reguła przypisująca im prawdopodobieństwa;
5. kalibracja wyniku liczbowego do wielkości i jednostki;
6. przynajmniej jedno przewidywanie możliwe do obalenia, nieużyte do wyznaczenia parametrów.

**[Udowodnione]** Powyższa lista nie zakłada kwantowości. W teorii klasycznej stanem może być miara na przestrzeni fazowej, instrumentem funkcja odczytu, a regułą prawdopodobieństwa push-forward tej miary. W teorii kwantowej naturalnymi odpowiednikami są macierze gęstości oraz instrumenty całkowicie dodatnie/POVM. W obu przypadkach sam generator nie określa przygotowania ani odczytu.

**[Silne dowody]** Minimalnym testem przyszłej fizycznej realizacji FIN powinno być: jedna jawna kalibracja oraz co najmniej dwa wcześniej niewykorzystane, bezwymiarowe ilorazy przewidziane poza próbą kalibracyjną. Dopasowanie samej skali do jednego wyniku nie jest przewidywaniem.

---

# 3. Faza 1 — niezależna rekonstrukcja minimalnego rdzenia

## 3.1. Dokładna algebra na \(C_{12}\)

Niech

\[
\mathcal H=\mathbb C^{12},\qquad
L=2I-T-T^*,
\]

gdzie \(T\) jest przesunięciem cyklicznym. W bazie Fouriera

\[
L e_k=\mu_k e_k,\qquad
\mu_k=2-2\cos\frac{2\pi k}{12}.
\]

Istnieje siedem różnych wartości \(\mu_0,\ldots,\mu_6\), ponieważ
\(\mu_k=\mu_{12-k}\).

**Twierdzenie 3.1 — dokładna algebra filtrów grafowych [Udowodnione].**  
Każda rzeczywista symetryczna macierz cyrkulantna, której współczynniki zależą wyłącznie od odległości cyklicznej, ma jednoznaczną reprezentację

\[
K=p(L),\qquad \deg p\le 6.
\]

**Dowód.** Radialna macierz na \(C_{12}\) ma siedem niezależnych wartości odpowiadających odległościom \(0,\ldots,6\), więc przestrzeń takich macierzy ma wymiar siedem. Jest diagonalna w bazie Fouriera, a jej symbol przyjmuje siedem wartości \(\kappa_0,\ldots,\kappa_6\). Węzły \(\mu_0,\ldots,\mu_6\) są parami różne. Interpolacja Lagrange’a daje jedyny wielomian \(p\) stopnia najwyżej sześć spełniający \(p(\mu_j)=\kappa_j\). Z rachunku funkcyjnego wynika \(K=p(L)\). \(\square\)

**[Udowodnione]** Dla zamrożonych profili strict i legacy obliczenie numeryczne dało względne reszty odpowiednio

\[
2.53\times10^{-14},\qquad 1.28\times10^{-13},
\]

czyli poziom błędu zmiennoprzecinkowego. Kod reprodukcyjny znajduje się w pliku
[fin_missing_principle_checks.py](./fin_missing_principle_checks.py).

**Wniosek 3.2 [Udowodnione].**  
Oba jądra, ich wspólne projekcje spektralne i wszystkie funkcje borelowskie należą do przemiennej algebry

\[
\mathcal A_{\rm rad}=C^*(L)\cong\mathbb C^7.
\]

**Wniosek falsyfikacyjny 3.3 [Udowodnione].**  
Fakt \(K_{\rm strict}=p_s(L)\), \(K_{\rm legacy}=p_l(L)\) nie wyróżnia tych dwóch jąder spośród wszystkich radialnych jąder \(C_{12}\). Jest dokładnym twierdzeniem reprezentacyjnym, ale sam nie jest prawem fizycznym ani dowodem unikalności.

## 3.2. Najmniejszy obiekt generujący wymagane konstrukcje

**Definicja 3.4 [Udowodnione jako definicja minimalna w podanej kategorii].**  
Dwuoznaczonym spektralnym układem dynamicznym nazywam

\[
\Sigma=(\mathcal H,L,F;\mathcal M,V),
\]

gdzie:

- \(\mathcal H\) jest przestrzenią Hilberta;
- \(L=L^*\) jest samosprzężonym operatorem bazowym;
- \(F:\{\mathrm{s},\mathrm{l}\}\times\sigma(L)\to\mathbb R\) jest parą oznaczonych filtrów;
- \(\mathcal M\) jest dopuszczalną rozmaitością/przestrzenią operatorów;
- \(V\) jest polem wektorowym określającym prawo adaptacji.

Konstrukcje są wtedy

\[
\begin{aligned}
K_s&=\int_{\sigma(L)}F(s,\lambda)\,dE_L(\lambda),\\
K_l&=\int_{\sigma(L)}F(l,\lambda)\,dE_L(\lambda),\\
G_L(z)&=\int_{\sigma(L)}\frac{1}{\lambda-z}\,dE_L(\lambda),\\
\dot K&=V(K).
\end{aligned}
\]

Miara spektralna \(E_L\) nie jest dodatkowym aksjomatem: daje ją twierdzenie spektralne.

### Test usuwania aksjomatów

| Składnik | Konieczność | Co dzieje się po usunięciu | Status |
|---|---:|---|---|
| Struktura Hilberta i \(L=L^*\) | 100/100 | znika kanoniczna PVM i zwykły samosprzężony rachunek spektralny | **Udowodnione** |
| Oznaczenie \(F(s,\cdot)\) | 95/100 | z algebry nie da się wskazać, który z nieskończenie wielu filtrów jest strict | **Udowodnione** |
| Oznaczenie \(F(l,\cdot)\) | 95/100 | analogicznie nie da się wskazać legacy | **Udowodnione** |
| Pole \(V\) | 100/100 dla prawa adaptacji | ten sam operator dopuszcza \(V=0\), przepływ gradientowy, unitarny, cieplny i wiele innych | **Udowodnione** |
| Skończony wymiar \(12\) | 15/100 abstrakcyjnie | nie jest potrzebny dla twierdzenia spektralnego; jest potrzebny dla konkretnego \(C_{12}\) i dokładnego rachunku siedmiowymiarowego | **Udowodnione** |
| Funkcjonał akcji i metryka gradientu | 0/100 dla ogólnego \(V\) | są potrzebne tylko, gdy żąda się dynamiki wariacyjnej, nie dowolnej adaptacji | **Udowodnione** |
| Interpretacja informacyjna | 0/100 dla rdzenia | nie zmienia operatora ani jego rachunku funkcjonalnego | **Udowodnione** |

**Twierdzenie 3.5 — brak redukcji do „gołego” operatora [Udowodnione].**  
Nie istnieje konstrukcja zależna wyłącznie od \(L\), która logicznie wymusza konkretną parę \((K_s,K_l)\) i konkretne \(V\). Dla tego samego \(L\) można wybrać dowolne dwie funkcje \(f_s,f_l\) na \(\sigma(L)\) oraz dowolne lokalnie Lipschitzowskie pole \(V\). Wszystkie wybory mają ten sam operator bazowy, a różne wyniki.

**[Udowodnione]** Można sztucznie zakodować \(F\) i \(V\) w jednym większym operatorze blokowym, lecz jest to zmiana notacji, nie redukcja zawartości aksjomatycznej. Dekodowanie bloków wymaga tych samych znaczników.

## 3.3. Co zostało z hipotezy „jednego obiektu”

**[Udowodnione]** Słaba wersja hipotezy przetrwała:

> strict, legacy, miara spektralna i rodzina Greena są różnymi funkcjami jednego samosprzężonego \(L\).

**[Udowodnione]** Silna wersja została sfalsyfikowana:

> samo istnienie \(L\) wybiera dokładne profile, adaptację, selektor i fizykę.

Przeciwprzykładem jest dowolna inna para filtrów i dowolne inne pole \(V\) na tym samym \(L\).

**[Silne dowody]** Najgłębszym czysto matematycznym opisem obecnego rdzenia jest zatem „oznaczony przemienny rachunek funkcyjny na skończonym operatorze grafowym”, a nie kompletna teoria pola, system wnioskujący ani fizyczna teoria informacji.

---

# 4. Faza 2 — czy selektor może wyłonić się z głębszej matematyki?

## 4.1. Twierdzenie o parzystej algebrze

Niech

\[
(R\psi)(x)=\psi(-x)
\]

będzie inwersją \(C_{12}\). Radialność daje

\[
RLR=L,\qquad RK_sR=K_s,\qquad RK_lR=K_l.
\]

Naturalny prąd orientacyjny

\[
J_{\rm cyc}=\frac{T-T^*}{2i}
\]

spełnia natomiast

\[
RJ_{\rm cyc}R=-J_{\rm cyc}.
\]

**Twierdzenie 4.1 — parzysta algebra nie tworzy nieparzystego źródła [Udowodnione].**  
Jeżeli

\[
\mathcal A_{\rm even}=C^*(I,L,K_s,K_l),
\]

to każdy \(A\in\mathcal A_{\rm even}\) spełnia \(RAR=A\), a zatem

\[
\mathcal A_{\rm even}\cap\{B:RBR=-B\}=\{0\}.
\]

Dotyczy to wielomianów, resolwent, semigrup cieplnych, projekcji spektralnych, funkcji zeta, determinantów oraz stanów postaci \(\rho=f(L,K_s,K_l)\).

**Dowód.** Generatory komutują z \(R\), więc komutuje z nim domknięta algebra przez nie generowana. Gdy \(B\) jest jednocześnie parzyste i nieparzyste, \(B=RBR=-B\), czyli \(B=0\). \(\square\)

## 4.2. Główne no-go selektora

Niech

\[
\mathcal O=\{+\omega,-\omega\}
\]

będzie torsorem orientacji, na którym \(R\omega=-\omega\).

**Twierdzenie 4.2 — brak naturalnego selektora [Udowodnione].**  
Dla danych \(D\) ustalonych przez inwersję \(R\) nie istnieje \(C_2\)-ekwiwariantna mapa

\[
S:D\longrightarrow\mathcal O.
\]

**Dowód.**

\[
S(D)=S(RD)=R\,S(D)=-S(D),
\]

lecz \(\mathcal O\) nie ma punktu stałego. \(\square\)

**[Udowodnione]** Dowód wymaga tylko rzeczywistej symetrii \(x\mapsto-x\), nie pełnej grupy \(\operatorname{Aut}(\mathbb Z_{12})\). Żaden algorytm diagonalizacji nie obchodzi przeszkody: równość \(\lambda_k=\lambda_{12-k}\) oznacza, że projekcja spektralna widzi przestrzeń \(\operatorname{span}\{e_k,e_{-k}\}\), ale nie wyróżnia \(+k\).

## 4.3. Bifurkacja i spontaniczne łamanie symetrii

**Twierdzenie 4.3 — symetryczny stan pozostaje symetryczny [Udowodnione].**  
Jeśli \(F\) jest lokalnie Lipschitzowskim polem \(C_2\)-ekwiwariantnym, \(F(Rx)=RF(x)\), i \(Rx_0=x_0\), to rozwiązanie

\[
\dot x=F(x),\qquad x(0)=x_0
\]

spełnia \(Rx(t)=x(t)\). Wynika to z jednoznaczności zagadnienia początkowego.

**[Udowodnione]** Potencjał

\[
V_a(m)=\frac14m^4-\frac a2m^2
\]

dla \(a>0\) ma dwa minima \(\pm\sqrt a\). Teoria bifurkacji wymusza parę gałęzi, nie jej element. Człon \(-hm\) wybiera znak, ale \(h\) jest nowym nieparzystym źródłem. Symetryczny szum może losowo wybrać gałąź w pojedynczej realizacji, lecz bez dodatkowego biasu rozkład znaków pozostaje symetryczny.

**[Udowodnione]** W pojedynczym układzie \(C_{12}\) nie istnieje limit termodynamiczny \(N\to\infty\), w którym wybór czystej fazy można zdefiniować kolejnością granic \(N\to\infty\), a następnie \(h\to0^\pm\). Taki mechanizm wymaga nowej rodziny układów.

## 4.4. Audyt kandydatów na ukryty selektor

| Mechanizm | Co rzeczywiście daje | Dlaczego nie wybiera znaku z obecnych danych | Status |
|---|---|---|---|
| Ekwariantna bifurkacja | istnienie par gałęzi o mniejszej symetrii | gałęzie są wymieniane przez \(R\) | **Udowodnione** |
| Spectral flow | całkowity indeks zorientowanej ścieżki | odwrócenie ścieżki zmienia znak; orientacja ścieżki jest wejściem | **Udowodnione** |
| Indeks i \(\eta\)-invariant | klasy stabilne po zadaniu operatora Diraca/gradingu | statyczne dane parzyste nie generują nieparzystego znaku | **Udowodnione w obecnym modelu** |
| \(K\)-teoria | klasy projekcji/faz po zadaniu szczeliny i symetrii | dla \(\mathbb C^7\): \(K_0\simeq\mathbb Z^7\), \(K_1=0\) | **Udowodnione** |
| Topologiczne \(H^1(S^1)\) | dowodzi istnienia dwóch orientacji | inwersja wysyła generator \(n\mapsto-n\) | **Udowodnione** |
| Morse/Conley | klasyfikuje punkty krytyczne i zbiory niezmiennicze | partnerzy \(x,Rx\) mają ten sam indeks | **Udowodnione** |
| Tomita–Takesaki | modularny przepływ pary algebra–stan | algebra przemienna daje przepływ trywialny; stan nieparzysty byłby nowym datum | **Udowodnione** |
| Crossed product | reprezentuje sektor znakowy | istnienie sektora nie wybiera obsadzenia ani wektora | **Udowodnione** |
| RG | przepływ sprzężeń | ekwiwariantny RG zachowuje podprzestrzeń parzystą | **Udowodnione** |
| Anomalia chiralna | może naruszyć klasyczną symetrię | wymaga fermionów, miary, regulatora i niezerowego indeksu | **Silne dowody: brak tych danych** |
| Sheaf/groupoid/higher groupoid | poprawnie koduje torsor i dane sklejenia | sekcja globalna nadal jest dodatkowym wyborem | **Udowodnione** |
| Derived geometry/HoTT | zachowuje symetrię locus krytycznego/ilorazu | formalizm nie tworzy punktu stałego | **Silne dowody** |
| Losowe macierze/free probability | rozkład statystyczny w zadanym ensemble/limicie | ensemble i orientacyjny bias są wejściami | **Silne dowody** |

**[Udowodnione]** Spectral flow spełnia

\[
\operatorname{sf}(D_{1-t})=-\operatorname{sf}(D_t),
\]

więc transportuje orientację parametru, ale jej nie tworzy. \(K\)-teoria faz topologicznych wymaga ponadto szczeliny, poziomu Fermiego i klasy symetrii. Są to struktury dodatkowe, a nie konsekwencje jednego radialnego operatora.

**[Udowodnione]** Dla przemiennej algebry \(\mathbb C^7\) modularna grupa każdego wiernego stanu jest trywialna. Dla rozszerzenia \(M_{12}(\mathbb C)\) i stanu \(\rho\)

\[
\sigma_t^\rho(A)=\rho^{it}A\rho^{-it},
\]

ale \(\rho\) jest nowym datum. Jeśli \(\rho=f(K_s,K_l,L)\), to \([R,\rho]=0\), więc modularność nadal nie łamie inwersji.

## 4.5. Minimalne rozwiązania problemu selektora

**[Udowodnione]** Istnieją dokładnie trzy logicznie odmienne rozwiązania:

1. **Iloraz gauge:** jeśli wszystkie obserwable spełniają \(O(+\omega)=O(-\omega)\), należy utożsamić dwa znaki. Selektor nie jest wtedy brakującą fizyką.
2. **Selekcja statystyczna:** teoria może przewidywać rozkład na dwóch gałęziach, np. \(1/2,1/2\), bez deterministycznego znaku.
3. **Źródło relacyjne:** jeśli znak jest obserwowalny, trzeba dodać stan, brzeg lub niezerowy obiekt \(J\) spełniający \(RJR=-J\), wraz z polaryzacją sprzężenia.

**[Udowodnione]** Dla samego kierunku minimalny deterministyczny dodatek to jeden punkt dwuelementowego torsora, informacyjnie jeden bit:

\[
\mathfrak o\in
\bigl(H^1(C_{12};\mathbb R)\setminus\{0\}\bigr)/\mathbb R_{>0}.
\]

Jeżeli trzeba również wybrać absolutny początek, minimalny wybór jest punktem torsora ramek

\[
\operatorname{Fr}(C_{12})=\mathbb Z_{12}\times\{\pm1\},
\]

czyli jednym z 24 elementów. To nie jest wynik funkcji spektralnej.

---

# 5. Faza 3 — czy matematyka może wygenerować wymiary fizyczne?

## 5.1. Twierdzenie naturalności jednostek

Wybierzmy jako bazowe wymiary długość \(L\), czas \(T\) i działanie \(A\). Grupa zmian jednostek

\[
\mathcal U=(\mathbb R_{>0})^3
\]

działa na wielkości o ciężarze \(w=(w_L,w_T,w_A)\) przez

\[
q\longmapsto u_L^{w_L}u_T^{w_T}u_A^{w_A}q.
\]

W szczególności:

\[
\begin{array}{c|c}
\text{wielkość} & (w_L,w_T,w_A)\\ \hline
\ell &(1,0,0)\\
t &(0,1,0)\\
\mathcal S &(0,0,1)\\
E &(0,-1,1)\\
m &(-2,1,1)
\end{array}
\]

**Twierdzenie 5.1 — brak niezerowej skali z danych bezwymiarowych [Udowodnione].**  
Niech \(X\) będzie zbiorem danych bezwymiarowych, na którym \(\mathcal U\) działa trywialnie. Jeżeli \(Q:X\to\mathbb R\) jest ekwiwariantną wielkością o niezerowym ciężarze \(w\), to \(Q=0\).

**Dowód.** Dla każdego \(x\in X\) i \(u\in\mathcal U\),

\[
Q(x)=Q(ux)=u^wQ(x).
\]

Ponieważ \(w\ne0\), można wybrać \(u\) z \(u^w\ne1\); stąd \(Q(x)=0\). \(\square\)

**Wniosek 5.2 [Udowodnione].**  
Bezwymiarowy operator może przewidywać bezwymiarowe ilorazy, wykładniki i porządki, ale nie może naturalnie podać dodatniej długości, czasu, działania, masy ani energii w jednostkach laboratoryjnych.

**[Udowodnione]** Ustawienie \(\hbar=c=k_B=1\) nie usuwa problemu. Jest wyborem jednostek, który identyfikuje wymiary, lecz powrót do sekund, metrów, dżuli i kelwinów wymaga kalibracji.

## 5.2. Najkrótsza jawna kalibracja

**[Udowodnione]** Po podaniu trzech dodatnich skal

\[
\ell_*,\qquad t_*,\qquad \mathcal A_*,
\]

bezwymiarowe liczby mogą zostać przeliczone:

\[
\ell_{\rm phys}=\ell_*d,\qquad
t_{\rm phys}=t_*\tau,\qquad
S_{\rm phys}=\mathcal A_*s,
\]

\[
E_*=\frac{\mathcal A_*}{t_*},\qquad
m_*=\frac{\mathcal A_*t_*}{\ell_*^2}.
\]

Jeśli teoria wymaga relatywistycznej relacji, \(c=\ell_*/t_*\) albo niezależny pomiar \(c\) jest dalszym prawem kalibracyjnym. W teorii kwantowej zwykle ustala się \(\mathcal A_*=\hbar\), lecz \(\hbar\) jest wtedy wejściem empirycznym.

## 5.3. Heat kernel i wymiar spektralny

Dla skończonego dodatniego \(L\)

\[
Z(t)=\operatorname{Tr}e^{-tL}
=\sum_{j=1}^{N}e^{-t\lambda_j},
\qquad
d_s(t)=-2\frac{d\log Z(t)}{d\log t}.
\]

**Twierdzenie 5.3 — brak asymptotycznego wymiaru z pojedynczego grafu skończonego [Udowodnione].**  
Dla ustalonego skończonego grafu \(Z(t)\to N\) przy \(t\to0\), więc \(d_s(t)\to0\). Dla grafu spójnego \(Z(t)\to1\) przy \(t\to\infty\), więc również \(d_s(t)\to0\).

**[Udowodnione]** Dla \(C_{12}\) obliczenie dało maksimum pośrednie

\[
\max_t d_s(t)=1.217782\quad\text{przy}\quad t\approx0.851138,
\]

ale wartości graniczne są numerycznie zerowe. Pośredni garb nie jest asymptotycznym wymiarem continuum. Stabilny wymiar wymaga rodziny \(L_N\), skali oczka i kontrolowanej granicy.

## 5.4. Zeta, spectral action i dylaton

**[Udowodnione]** Dla skończonego dodatniego operatora

\[
\zeta_D(s)=\sum_j\lambda_j^{-s}
\]

spełnia

\[
\zeta_{cD}(s)=c^{-s}\zeta_D(s).
\]

Analityczna kontynuacja albo wyznacznik zeta nie wybiera \(c\); jedynie śledzi jego ciężar. W zastosowaniach z wymiarowym operatorem regularizacja często wymaga skali renormalizacji \(\mu\).

**[Udowodnione]** Zasada działania spektralnego ma postać

\[
\operatorname{Tr}f(D^2/\Lambda^2).
\]

Funkcja odcięcia \(f\) i skala \(\Lambda\) są dodatkowymi danymi. Dodanie dylatonu może zastąpić stałą skalę polem, ale wtedy pole, jego akcja, stan i wybór wartości oczekiwanej są nową strukturą.

## 5.5. Modular flow i czas termiczny

**[Udowodnione]** Hipoteza czasu termicznego daje kanoniczny przepływ względem pary \((\mathcal M,\rho)\), nie względem samej algebry. Dla \(\mathcal M=M_n(\mathbb C)\)

\[
\sigma_\tau^\rho(A)=\rho^{i\tau}A\rho^{-i\tau}.
\]

Jeżeli \(\rho=e^{-\beta H}/Z\), to

\[
\sigma_\tau^\rho(A)
=e^{-i\beta H\tau}Ae^{i\beta H\tau}.
\]

Porównanie z ewolucją Heisenberga pokazuje, że czas fizyczny jest proporcjonalny do \(\beta\hbar\tau\), z zależnością znaku od konwencji. Stan \(\rho\), \(\beta\hbar\) i zegar nie wynikają z samego \(K\).

**[Udowodnione]** W przemiennej algebrze \(C^*(K)\) przepływ modularny jest trywialny. Zatem Tomita–Takesaki nie dostarcza ukrytego czasu obecnemu przemiennemu rdzeniowi.

## 5.6. RG i transmutacja wymiarowa

**[Udowodnione]** Równanie RG

\[
\mu\frac{dg}{d\mu}=\beta(g)
\]

może wytworzyć niezmiennik

\[
\Lambda=\mu\exp\!\left(-\int^{g(\mu)}
\frac{dg}{\beta(g)}\right).
\]

Jest to ważna transmutacja bezwymiarowego sprzężenia w skalę **względem** \(\mu\) i warunku renormalizacji. Nie generuje absolutnej jednostki z bezwymiarowego obiektu; skala transformuje się wraz ze zmianą kalibracji. Ponadto wymaga rodziny teorii, regulatora i niezerowej funkcji beta, których pojedyncza macierz \(12\times12\) nie określa.

**[Silne dowody]** Mechanizm Colemana–Weinberga może dynamicznie wybrać niezerową wartość oczekiwaną, lecz wymaga pól kwantowych, oddziaływań, miary, renormalizacji i wyboru próżni. Jest kandydatem dla przyszłego rozszerzenia, nie twierdzeniem o obecnym operatorze.

## 5.7. Pozostałe drogi do geometrii i skali

| Droga | Co może wyłonić | Niezbędne dane, których nie daje sam operator | Status |
|---|---|---|---|
| Fisher information | lokalną metrykę na rodzinie rozkładów | rodzina statystyczna, parametr, eksperyment | **Udowodnione** |
| Information geometry | połączenia i krzywizny przestrzeni modeli | identyfikacja stanów i obserwacji | **Udowodnione** |
| Optimal transport | metrykę Wassersteina i przepływy gradientowe | przestrzeń bazowa, koszt/metryka, miary, mobilność | **Udowodnione** |
| Entropy production | strzałkę w zadanej dynamice otwartej | generator Markowa/Lindblada, stan odniesienia, czas | **Udowodnione** |
| Geometric flows | ewolucję metryki | początkowa geometria, równanie i parametr czasu | **Udowodnione** |
| Causal ordering | relację „przed–po” | gęstość/elementarna objętość dla skali metrycznej | **Silne dowody** |
| Causal sets | porządek i objętość przez zliczanie | sprinkling density lub fundamentalna skala dyskretności | **Silne dowody** |
| Tensor networks | geometrię korelacji w zadanej rodzinie stanów | tensory, bond dimensions, stan, sieć i granica | **Silne dowody** |
| Operator scaling limit | continuum i efektywny generator | sekwencja \(K_N\), \(a_N\), mapy porównawcze i zbieżność | **Udowodnione** |

**Twierdzenie 5.4 — pojedynczy \(K_{12}\) nie wyznacza continuum [Udowodnione].**  
Istnieją nieskończenie liczne sekwencje \(\{K_N\}\) przechodzące przez ten sam element \(K_{12}\), lecz z limitami lokalnymi, nielokalnymi, rozłącznymi albo bez limitu. Wartość jednego elementu ciągu nie determinuje całego diagramu refinement.

---

# 6. Faza 4 — most informacja → fizyka

## 6.1. Nieidentyfikowalność energii i temperatury

**Twierdzenie 6.1 — klasyczna nieidentyfikowalność Gibbsa [Udowodnione].**  
Niech \(p_i>0\), \(\sum_i p_i=1\). Dla dowolnych \(\theta>0\) i \(C\in\mathbb R\) połóżmy

\[
E_i^{(\theta,C)}=-\theta\log p_i+C.
\]

Wtedy

\[
p_i=
\frac{e^{-E_i^{(\theta,C)}/\theta}}
{\sum_j e^{-E_j^{(\theta,C)}/\theta}}.
\]

Ten sam rozkład, a więc ta sama informacja Shannona, jest zgodny z nieprzeliczalną rodziną skal energii, temperatur i zer energii.

**Twierdzenie 6.2 — kwantowa nieidentyfikowalność Gibbsa [Udowodnione].**  
Dla dowolnego wiernego stanu \(\rho\)

\[
H_{\theta,C}=-\theta\log\rho+CI
\]

daje \(\rho=e^{-H_{\theta,C}/\theta}/Z\) dla każdego \(\theta>0\). Hamiltonian modularny lub splątania \(-\log\rho\) jest bezwymiarowy; fizyczna energia wymaga \(\theta\).

**[Udowodnione]** Obliczenie kontrolne zachowało entropię jednego równomiernego bitu równą \(\ln2\), jednocześnie zmieniając średnią energię z \(0.5\) na \(50\) przez samą zmianę skali poziomów. Informacja nie identyfikuje energii.

## 6.2. Jaynes

**[Udowodnione]** Maksymalizacja

\[
H(p)=-\sum_i p_i\ln p_i
\]

przy ograniczeniach normalizacji i zadanej średniej energii daje

\[
p_i=Z^{-1}e^{-\beta E_i}.
\]

Metoda Jaynesa wyprowadza postać rozkładu po zadaniu poziomów \(E_i\), miary bazowej i ograniczeń. Nie wyprowadza z informacji tego, które wielkości są energią ani jaka jest ich fizyczna skala.

## 6.3. Dokładny most Landauera

Niech pamięć \(S\) początkowo nie koreluje z rezerwuarem

\[
\tau_R=\frac{e^{-\beta H_R}}{Z_R},
\]

a całość przechodzi ewolucję unitarną. Niech

\[
\Delta s=S(\rho_S)-S(\rho'_S)
\]

będzie spadkiem bezwymiarowej entropii pamięci, a
\(Q=\operatorname{Tr}H_R(\rho'_R-\tau_R)\) ciepłem oddanym do rezerwuaru.

**Twierdzenie 6.3 — skończeniewymiarowa tożsamość Landauera [Udowodnione warunkowo].**

\[
\beta Q
=\Delta s+I(S{:}R)_{\rho'}
+D(\rho'_R\Vert\tau_R)
\ge\Delta s.
\]

**Dowód.** Unitarność i początkowy iloczyn dają
\(\Delta S_R=\Delta s+I(S{:}R)_{\rho'}\). Z definicji stanu Gibbsa

\[
D(\rho'_R\Vert\tau_R)=\beta Q-\Delta S_R.
\]

Po złożeniu otrzymujemy tożsamość i nierówność. \(\square\)

**[Udowodnione warunkowo]** Dla wymazania jednego bitu

\[
Q\ge k_BT\ln2.
\]

Jeśli liczba \(4\ln2\) opisuje cztery wymazywane bity, to

\[
Q_{\min}=4k_BT\ln2.
\]

Nie jest to unikalna predykcja FIN: ten sam bound obowiązuje dowolną fizyczną pamięć z takim samym spadkiem entropii. Potrzebne są \(H_R\), kąpiel, \(T\), proces wymazania i identyfikacja logicznych stanów.

## 6.4. Względna entropia, KMS i teorie zasobów

**Twierdzenie 6.4 [Udowodnione warunkowo].**  
Dla \(\tau_\beta=e^{-\beta H}/Z\)

\[
D(\rho\Vert\tau_\beta)
=\beta\bigl(F_\beta(\rho)-F_\beta(\tau_\beta)\bigr),
\]

gdzie

\[
F_\beta(\rho)=\operatorname{Tr}(\rho H)
-\beta^{-1}S(\rho).
\]

Względna entropia staje się fizyczną różnicą energii swobodnej dopiero względem zadanego \(H\) i \(\beta\).

**[Udowodnione]** Warunek KMS definiuje równowagę względem zadanej algebry i zadanej dynamiki automorfizmów. Nie wytwarza generatora ani temperatury z samej algebry.

**[Udowodnione]** Termodynamiczne teorie zasobów zakładają Hamiltonian, temperaturę, stany termiczne i klasę dozwolonych operacji. Ich monotony są potężnymi twierdzeniami po tych założeniach, nie zasadą generującą założenia.

## 6.5. Werdykt informacyjny

**[Udowodnione]** Shannon i von Neumann mierzą bezwymiarową niepewność. Fizyczna entropia ma postać

\[
S_{\rm th}=k_B S_{\rm info}
\]

dopiero po fizycznej identyfikacji stanów, ensemble’u i kalibracji temperatury/energii. Stała \(k_B\) jest przelicznikiem jednostek i skalą empiryczną, nie wynikiem rachunku spektralnego.

**[Udowodnione]** Najkrótszy prawidłowy most brzmi:

\[
\text{informacja}
+\text{Hamiltonian}
+\text{stan/kąpiel}
+\text{temperatura}
+\text{dozwolony proces}
\Longrightarrow
\text{entropia i koszt fizyczny}.
\]

Usunięcie któregokolwiek z fizycznych wejść przywraca nieidentyfikowalność.

