# DKO-001–DKO-005 — wyjaśnienie po polsku

**O co chodzi?** Po wcześniejszym odkryciu, że różne operatory mogą generować te same zbiorowe fluktuacje krytyczne, kampania badała, czy *czas* przełączeń pozwala je rozróżnić bez mierzenia bezwzględnej amplitudy sygnału.

**Nowe narzędzie:** mierzymy średnią z wcześniej wybranej połowy elementów i rozdzielamy wynik na wspólną zmianę całej populacji F oraz dodatkową fluktuację bloku U. Wyprowadzono z prawdziwego generatora FIN (a nie założono) dwie relaksacje: F zanika w tempie γ, U w tempie 1. Powstają dwa różne składniki szumu granicznego. Dla jednego kanału ich suma pamięta przeszłość i zwykle nie zachowuje się jak prosty pojedynczy proces Markowa.

**Najważniejsza prognoza:** porównujemy tylko znaki sygnału w chwili 0 i 1. W granicznym opisie dla A7 znaki są jednakowe z prawdopodobieństwem ~78,6077%, a w operatorze porównawczym ~78,8475%. Drobna różnica jest ściśle dodatnia na mocy wzoru analitycznego, lecz nie jest jeszcze pomiarem fizycznym. Internalny tryb kontrolny r1 pozwala usunąć wspólną niewiadomą jednostkę czasu z porównania; nie usuwa błędów czujników i offsetów.

**Koszt obserwacji:** dla gaussowskiego procesu granicznego, jeśli obserwujemy udział α populacji M razy, siła rozróżnienia rośnie zasadniczo jak Mα². Przy Mα²→0 rozróżnienie jest niemożliwe; przy Mα²→∞ test znakowy potrafi rozróżniać z ryzykiem dążącym do zera. Nie oznacza to, że w prawdziwym skończonym układzie wystarczy dowolnie mały blok i bardzo długi czas: jednoczesnej granicy jeszcze nie udowodniono.

**Status:** warunkowa teoria matematyczna istniejącego modelu FIN, 32 niezależne kontrole algebry i skończonego generatora. Nadal brak dowodu, że przyroda realizuje A7, oraz danych dotyczących fizycznego czasu, aparatury i kalibracji.
