# FIN ST2352–ST2366 — audyt graniastosłupa prawidłowego trójkątnego

## Werdykt skrócony

Graniastosłup trójkątny jest **bardzo dobrym lokalnym obiektem dla refinacji FIN**, ponieważ jest dokładnie komórką produktową

\[
\Delta^2\times I.
\]

Łączy trójpunktową ścianę na jednej warstwie z odpowiadającą jej ścianą na następnej warstwie. Trzy ściany boczne przenoszą trzy krawędzie trójkąta.

Nie jest jednak brakującym źródłem fizyki. Nie wybiera:

- trójkąta bazowego,
- kierunku między warstwami,
- nowej informacji trójpunktowej,
- stosunku wag Hodge’a,
- skali długości ani czasu,
- obserwatora lub aparatury.

Najlepsza rola naukowa brzmi:

> **kanoniczna warunkowa komórka lokalna transportująca już istniejącą strukturę trójkątną pomiędzy dwiema wybranymi warstwami.**

## 1. Definicja geometryczna

Prawidłowy prosty graniastosłup trójkątny ma dwie równoległe podstawy będące trójkątami równobocznymi. Dla boku podstawy \(a\) i wysokości \(h\) można wybrać wierzchołki

\[
\begin{aligned}
v_0&=(0,0,-h/2),&v_3&=(0,0,h/2),\\
v_1&=(a,0,-h/2),&v_4&=(a,0,h/2),\\
v_2&=(a/2,\sqrt3a/2,-h/2),&v_5&=(a/2,\sqrt3a/2,h/2).
\end{aligned}
\]

Objętość i pole powierzchni wynoszą

\[
V=\frac{\sqrt3}{4}a^2h,
\qquad
S=\frac{\sqrt3}{2}a^2+3ah.
\]

Określenie „prawidłowy” nie musi oznaczać \(h=a\). Warunek \(h=a\) daje szczególny graniastosłup jednostajny: wszystkie dziewięć krawędzi ma tę samą długość, a ściany boczne są kwadratami.

## 2. Dokładna struktura komórkowa

Graniastosłup ma:

| wymiar | komórki | liczba |
|---:|---|---:|
| 0 | wierzchołki | 6 |
| 1 | krawędzie | 9 |
| 2 | dwie podstawy trójkątne i trzy ściany boczne | 5 |
| 3 | wnętrze graniastosłupa | 1 |

Jest to dokładnie rozkład produktu wypełnionego trójkąta i odcinka:

\[
\begin{aligned}
C_0&: 3\cdot2=6,\\
C_1&: 3\cdot2+3\cdot1=9,\\
C_2&: 1\cdot2+3\cdot1=5,\\
C_3&: 1\cdot1=1.
\end{aligned}
\]

Dla orientacji przyjętej w certyfikacie druga macierz brzegowa ma postać

\[
d_1=
\begin{pmatrix}
-1&-1& 1& 0& 0& 0& 0& 0& 0\\
 0& 0& 0& 1& 1&-1& 0& 0& 0\\
 1& 0& 0&-1& 0& 0&-1& 1& 0\\
 0& 1& 0& 0&-1& 0& 0&-1& 1\\
 0& 0&-1& 0& 0& 1& 1& 0&-1
\end{pmatrix}.
\]

Dokładnie:

\[
d_1d_0=0,
\qquad
d_2d_1=0,
\]

oraz

\[
\operatorname{rank}d_0=5,
\qquad
\operatorname{rank}d_1=4,
\qquad
\operatorname{rank}d_2=1.
\]

Sama powierzchnia jest sferą \(S^2\), więc

\[
(b_0,b_1,b_2)=(1,0,1).
\]

Po dodaniu komórki objętościowej otrzymujemy kulę \(B^3\):

\[
(b_0,b_1,b_2,b_3)=(1,0,0,0).
\]

## 3. Widma

Graf graniastosłupa jest iloczynem kartezjańskim

\[
C_3\square K_2.
\]

Widmo jego laplasjanu grafowego wynosi

\[
\operatorname{spec}L_0=\{0,2,3,3,5,5\}.
\]

Dla nieważonego komórkowego laplasjanu jednoforem

\[
L_1=d_0d_0^*+d_1^*d_1
\]

otrzymujemy dokładnie

\[
\operatorname{spec}L_1=
\{2,3,3,3,5,5,5,5,5\}.
\]

Nie ma jednoforem harmonicznych, zgodnie z \(H^1=0\).

## 4. Dwie niewyznaczone wagi Hodge’a

Symetria graniastosłupa zachowuje dwa typy ścian:

1. dwie ściany trójkątne,
2. trzy ściany kwadratowe/prostokątne.

Niech \(t>0\) będzie wagą ścian trójkątnych, a \(q>0\) wagą ścian bocznych. Wtedy sektor gradientowy ma stałe widmo

\[
\{2,3,3,5,5\},
\]

a sektor koexact ma dokładne wartości

\[
\boxed{
3t,
\quad
5q\ \text{(dwukrotnie)},
\quad
3t+2q.}
\]

Graniastosłup nie wybiera stosunku \(t/q\). Nawet geometryczny warunek \(h=a\) nie wystarcza automatycznie, ponieważ dyskretna gwiazda Hodge’a zależy od przyjętej konstrukcji dualnej, nie tylko od pól ścian.

## 5. Pełny skan osadzeń w strict FIN

Na 12 oznaczonych wierzchołkach istnieje dokładnie

\[
55\,440
\]

różnych podgrafów będących graniastosłupami trójkątnymi.

Przebadano wszystkie osadzenia, stosując dwa naturalne wyniki jakości:

- sumę dziewięciu wag strict \(W_{ij}\),
- iloczyn dziewięciu wag, obliczany przez sumę logarytmów.

| kryterium | najlepszy wynik | degeneracja | luka do następnego poziomu |
|---|---:|---:|---:|
| suma wag | 2.8495753042468612 | 12 | 0.021501500150249875 |
| logarytm iloczynu | -13.189306060622824 | 12 | 0.002470968459169498 |

W obu przypadkach 12 maksimów tworzy jedną pełną orbitę \(D_{12}\), ale oba kryteria wybierają **różne** graniastosłupy.

Wynika z tego:

\[
\boxed{
\text{strict }W
\text{ nie wybiera kanonicznego graniastosłupa bez wyboru funkcjonału.}}
\]

Nawet po wybraniu funkcjonału pozostaje 12-krotna degeneracja translacyjno-odbiciowa.

## 6. Symetria i kierunek warstw

Grupa symetrii prawidłowego graniastosłupa trójkątnego ma rząd 12:

\[
D_{3h}\simeq D_3\times C_2.
\]

Zawiera transformację wymieniającą podstawę dolną z górną. Dlatego bryła sama nie rozstrzyga:

- która warstwa jest wcześniejsza,
- która jest późniejsza,
- jaki jest znak dyskretnej pochodnej między warstwami.

Uznanie jednej podstawy za „wejście”, a drugiej za „wyjście” jest dodatkowym zegarem/orientacją. Graniastosłup nie usuwa przeszkody `QW-2191`.

## 7. Co graniastosłup daje FIN

### Wyniki dodatnie

- Daje kanoniczną komórkę produktu dla **już wybranego** trójkąta.
- Dostarcza dwie kopie ściany trójkątnej i trzy konieczne ściany pionowe.
- Przy pierwotnej całkowitej incydencji usuwa dowolność znaków/współczynników operatora brzegowego.
- Zapewnia dokładne tożsamości \(d^2=0\).
- Może transportować \(\tau_{ijk}\) lub parzystość \(x_1x_2x_3\) między warstwami.
- Komórka objętościowa daje naturalne miejsce dla operatora \(d_2\) i dyskretnego prawa zamknięcia.

### Czego nie daje

- Nie generuje nowego połączonego kumulantu trzeciego rzędu.
- Nie wybiera trójkąta bazowego w pełnym grafie strict.
- Nie wybiera jednego z 55 440 osadzeń.
- Nie wybiera kierunku warstw.
- Nie wybiera wag \(t,q\).
- Nie ustala \(a,h\) ani ich jednostek.
- Nie tworzy stożka przyczynowego, teorii Maxwella ani aparatury.

## 8. Relacja z 12-stanowym nośnikiem

Pojedynczy graniastosłup trójkątny ma sześć wierzchołków, więc nie jest bezpośrednio całym nośnikiem \(Z_{12}\).

Możliwe konstrukcje 12-wierzchołkowe, takie jak dwa graniastosłupy albo graniastosłup sześciokątny, wymagają dodatkowego rozkładu. Nie jest on wyznaczony przez samą strukturę cykliczną FIN. W szczególności dwa rozłączne graniastosłupy nie mają topologii ani widma cyklu \(Z_{12}\).

Najbardziej naturalne użycie pozostaje lokalne:

\[
\text{wybrany trójkąt FIN}
\quad\longmapsto\quad
\text{jego graniastosłup transportowy }\Delta^2\times I.
\]

## 9. Konkluzja

Graniastosłup prawidłowy trójkątny jest matematycznie właściwą odpowiedzią na pytanie:

> Jak transportować jedną trójkątną komórkę informacyjną między dwiema warstwami refinacji, zachowując pełną strukturę brzegową?

Nie odpowiada natomiast na pytanie:

> Dlaczego FIN ma wybrać właśnie ten trójkąt, tę orientację, ten stosunek wag i tę fizyczną skalę?

Werdykt:

\[
\boxed{
\text{ważny lokalny obiekt konstrukcyjny — nie brakująca zasada fundamentalna.}}
\]

## 10. Rozszerzenie na pełny kompleks strict

Jeżeli komórkę \(\Delta^2\times I\) zastosujemy nie do jednego, lecz do wszystkich 220 trójkątów pełnego kompleksu flagowego strict, otrzymujemy produkt

\[
\Delta^{11}\times I.
\]

Liczba komórek w wymiarach \(0,\ldots,12\) wynosi dokładnie

\[
(24,144,506,1210,2079,2640,2508,1782,935,352,90,14,1).
\]

W szczególności:

- istnieje 220 nakładających się graniastosłupów trójkątnych,
- istnieje 66 pionowych kwadratów — po jednym dla każdej krawędzi bazowej,
- kwadraty tworzą 6 orbit odległościowych \(D_{12}\),
- całkowita pierwotna incydencja usuwa dowolność ich współczynników brzegowych.

Pełny produkt jest kontrakcyjny, lecz jego sam 2-szkielet ma

\[
\dim H^2=385,
\]

dopóki nie zostaną dołączone komórki wyższych wymiarów. Konstrukcja jest więc algebraicznie spójna, ale bardzo duża i nielokalna.

## 11. Niestabilność wyboru funkcjonału

Dla rodziny wyników

\[
S_\alpha(P)=\sum_{e\in P}W_e^\alpha
\]

dwie najlepsze orbity graniastosłupów zamieniają kolejność przy

\[
\boxed{\alpha_*=0.0055703936174915595.}
\]

W punkcie przejścia degeneracja wynosi 24, ponieważ remisują dwie pełne orbity \(D_{12}\). Po obu stronach pozostaje degeneracja 12. Jest to bezpośredni dowód, że niewielka zmiana niesourced funkcjonału może zmienić wybrany motyw.

## 12. Ostateczny gate ST2352–ST2441

| wymaganie | wynik |
|---|---|
| dokładna komórka \(\Delta^2\times I\) | PASS |
| pierwotna incydencja i \(d^2=0\) | PASS |
| warunkowe domknięcie ścian pionowych | PASS |
| strict selektor trójkąta bazowego | FAIL |
| kanoniczne globalne osadzenie | FAIL |
| kierunek warstw/czasu | FAIL |
| nowe źródło trójpunktowe | FAIL |
| unikalny stosunek Hodge’a | FAIL |
| fizyczna skala i continuum | FAIL |
| OA/aparatura/rekord | FAIL |

Końcowy wynik:

\[
\boxed{3/10\ \text{warunkowych wierszy matematycznych},\qquad 0/10\ \text{wierszy fizycznych}.}
\]

## Artefakty wykonawcze

- `fin_st2352_st2366_triangular_prism_audit.py`
- `test_fin_st2352_st2366_triangular_prism.py`
- `FIN_ST2352_ST2366_Results.json`
- `FIN_ST2352_ST2366_Summary.csv`
- pakiety `FIN_ST2352_*.json`–`FIN_ST2366_*.json`
- `fin_st2367_st2441_triangular_prism_followup.py`
- `test_fin_st2352_st2441_prism_cycle.py`
- zbiory wyników `FIN_ST2367_ST2381_Results.json`–`FIN_ST2427_ST2441_Results.json`
- `FIN_ST2427_ST2441_TriangularPrismGate.json`
