# Wyprowadzone sprzężenie i konkretna relacja szybkości FIN

Data: 2026-10-16. Status: **dokładna relacja warunkowa w przyjętym modelu;
niezależne źródło kinetyki pozostaje otwarte**.

## Wynik

Z rozkładu Gibbsa FIN i szczegółowej równowagi można wyprowadzić sprzężenie
dwóch zmian etykiet, które narzuca relację między szybkościami przejść.
Dla aktualnego A7 otrzymuję

\[
\boxed{\frac{\log R_3}{\log R_6}
\in[0.77490314788261,\;0.77490314788262].}
\]

Poniżej definiuję R dokładnie. Relacja obowiązuje przy tym samym N i g>0,
ale jej wartość nie zależy od tych parametrów, wspólnej jednostki czasu ani
wyboru heat-bath, Metropolisa czy Barkera. Wymaga natomiast zadanej postaci
Gibbsa, szczegółowej równowagi i właściwie określonych szybkości etykietowanych.

Nie jest to wyprowadzenie nowego mikroskopowego prawa. Statyczny odpowiednik
sprzężenia występuje już w RG-E02. Tutaj zapisano jego konsekwencję dla
dowolnej odwracalnej kinetyki, poprawkę dla liczności oraz certyfikowaną
relację liczbową dla A7.

## 1. Skąd wynika sprzężenie

Stan x opisuje N etykiet; n(x) oznacza ich liczności. Przyjmujemy istniejące

\[
\pi_N(x)\propto
\exp\!\left[\frac{g}{2N}n(x)^T A_7n(x)+\sum_i h_i(x_i)\right].
\]

Pola h_i mogą być niezerowe, ale nie mogą zmieniać się wraz z kontekstem
pozostałych etykiet. Szczegółowa równowaga daje dla rzeczywistych szybkości

\[
\frac{q(x,y)}{q(y,x)}=\frac{\pi_N(y)}{\pi_N(x)}.
\]

Niech pierwsza wybrana etykieta zmienia się z a na b. Drugą, odrębną,
zmieniamy z c na d; resztę z pozostawiamy taką samą. Definiuję

\[
r_{ab}(c;z)=
\frac{q((a,c,z),(b,c,z))}{q((b,c,z),(a,c,z))}.
\]

Odjęcie dwóch logarytmów i rozwinięcie formy kwadratowej daje

\[
\boxed{\mathcal C_{ab;cd}
=\log\frac{r_{ab}(d;z)}{r_{ab}(c;z)}
=\frac gN(A_{bd}-A_{ad}-A_{bc}+A_{ac}).}
\]

Wkład pozostałych etykiet, stałych pól i wyrazów własnych kasuje się.
To dokładny wynik dla skończonego N, bez granicy asymptotycznej i bez
dopasowania szybkości. W zapisie wektorowym wyprowadzone sprzężenie ma postać

\[
\mathcal B_N(u,v)=\frac gN u^T A_7v,
\qquad u=e_b-e_a,\quad v=e_d-e_c.
\]

Ponieważ A7=XX^T, sprzężenie jest symetryczne, a dla tej samej zmiany

\[
\mathcal B_N(u,u)=\frac gN\lVert X_b-X_a\rVert^2\ge0.
\]

Dotyczy to przestrzeni cech modelu. Nie nadaje etykietom roli punktów
fizycznej przestrzeni i nie tworzy nowej warstwy pod nadsolitonem.

## 2. Relacja, którą można obliczyć bez znajomości aktywności kinetycznej

Dla d=3 i d=6 wybieramy a=c=0 oraz b=d. Oznaczam

\[
R_d=\frac{r_{0d}(d;z)}{r_{0d}(0;z)}.
\]

Jest to iloraz dwóch stosunków szybkości: pierwszy przy drugiej etykiecie
równej d, drugi przy drugiej etykiecie równej 0. Każdy stosunek używa
przejścia i jego dokładnego odwrócenia. Wtedy

\[
\log R_d=\frac{2g}{N}(A_0-A_d).
\]

Z istniejącego rozkładu Fouriera A7 wynika

\[
A_0-A_3=\frac{\Lambda_3+\Lambda_5+\Lambda_6}{6},\qquad
A_0-A_6=\frac{\Lambda_3+\Lambda_5}{3}.
\]

Stąd

\[
\boxed{R_3=R_6^{\rho_{36}},\qquad
\rho_{36}=\frac12\left(1+\frac{\Lambda_6}{\Lambda_3+\Lambda_5}\right).}
\]

Dla zadanej w repozytorium matematycznej postaci kernela
`cos(0.18575 d+0.1625)/(1+d^(9/5))` i projekcji A7 przedział podany na
początku jest certyfikowany od danych wejściowych do wyniku. Obliczenie używa
liczb wymiernych, ograniczeń szeregu cosinusa i pierwiastków sprawdzanych
nierównościami całkowitoliczbowymi. Nie opiera się na traktowaniu liczby
zmiennoprzecinkowej jako przedziału.

Istniejący kontrmodel `FLAT_P7_TRACE`, o równych wartościach własnych na tym
samym siedmiowymiarowym nośniku, daje dokładnie `rho_36=3/4`. Jest więc
rzeczywista różnica między przewidywaniami modeli. Nie jest to jeszcze wynik
pomiaru ani argument, że natura wybiera A7.

Przy g=0 oba R są równe 1, więc iloraz ich logarytmów jest nieokreślony.
Blisko zera potrzebny byłby osobny budżet błędów pomiaru. Numerów 3 i 6 nie
należy interpretować jako odległości w jednostkach fizycznych.

## 3. Ważna poprawka dla procesu liczności

Jeśli korzystamy z liczności zamiast konkretnej etykiety, przejście a->b
może wykonać n_a kopii. Przejście odwrotne ma już n_b+1 możliwych kopii.
Właściwa wielkość to dlatego

\[
\log\frac{Q(n,n-e_a+e_b)}{Q(n-e_a+e_b,n)}
-\log\frac{n_a}{n_b+1}.
\]

Bez tej poprawki wynik jest zafałszowany. W sprawdzonym przykładzie dwóch
kolejnych zmian 0->3 surowy kontrast różni się od przewidywanego sprzężenia
dokładnie o `−log(3)`. Ponadto warunkowanie na pełnym stanie nie jest
automatycznie dostępne z instrumentu blokowego MM.

## 4. Co nadal pozostaje swobodne

Ta relacja ogranicza kombinacje szybkości, lecz nie wyznacza całego generatora.
Pokazuje to jawny kontrmodel: każdą szybkość heat-bath zmiany o wewnętrzny
dystans etykiet równy 3 mnożymy przez `1/(1+g^2)`; pozostałe zachowujemy.

Mnożnik jest taki sam w obu kierunkach, więc zachowuje rozkład Gibbsa,
szczegółową równowagę, symetrie i wszystkie powyższe R. Może działać przez
odrzucanie propozycji przy tym samym jednostkowym zegarze prób. Przy g=0
proces jest dokładnie wyjściowym heat-bath. Przy g=3 stosunek szybkości
wychodzących 0->3 oraz 0->6 zmienia się jednak dziesięciokrotnie.

To kontrprzykład analityczny, nie przyjęty nowy model FIN. Pokazuje dokładnie,
dlaczego otrzymana predykcja nie jest jeszcze niezależnym źródłem kinetyki.
Nadal brakuje zasady określającej symetryczne aktywności przejść.

Nie ma tu sprzeczności z wcześniejszym W=0: tam wykluczono dodatkowy
odwracalny przepływ w niezmienionym działaniu obsadzeń. Tutaj wyprowadzono
sprzężenie sił określających stosunki szybkości. W równowadze prądy netto
pozostają zerowe.

## 5. Weryfikacja i status

[Dowód](PROOF.md) zawiera pełne przesłanki, zakres i kontrprzykład.
[Certyfikat](CERTIFICATE.json) rozdziela dokładne ograniczenia od
[diagnostyk numerycznych](RESULTS.json). Przeszło 576 dokładnych przypadków
algebraicznych oraz osiem małych porównań szybkości w czterech kinetykach
i dwóch kontekstach. Te przypadki kontrolne nie zastępują dowodu ogólnego.

Nie zmieniono `AGENTS.md`, wcześniejszych pakietów ani podstawowego prawa FIN.
Wynik jest warunkową predykcją modelu. **K07 nadal nie został dostarczony;
KB-002/003 pozostają nieaktywne.** Pomiar wymagałby osobnego uprawnionego
interfejsu i jednoczesnych przedziałów błędu, zamiast bramki dokładnej równości.
