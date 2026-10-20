# Wyprowadzona granica dynamiki z hipotezy historii

Data: 2026-10-20. Źródło: niezmieniona, wcześniej upoważniona hipoteza
FIN-HYP-HIST-ME-001. Nie jest to strict FIN ani wyprowadzone prawo natury.

## Nowy wynik dodatni

**Z konkretnej rodziny źródłowej K+cI wyprowadziłem generator graniczny**,
zamiast wstawiać go jako równanie wejściowe:

    Q_xy=K_xy sqrt(pi_y/pi_x), x!=y,
    Q_xx=-sum_(y!=x)Q_xy.

K oznacza tu macierz legalnych kontynuacji konfiguracji, nie jądro A7
ani kernel strict/legacy; oddziaływanie A7 występuje w zadanej mierze pi.

Znana postać pierwiastkowych szybkości nie jest nowym odkryciem ogólnym.
Nowy krok lokalny to opłacona droga od zadanego skończonego źródła,
przez pominięcie jednostajnych współrzędnych, do generatora i prawa
trajektorii — z jawnymi błędami i bez dopasowywania go pod HB.

Przy stałym skończonym nośniku i dodatniej pi udowodniono błąd skalowania
normalizującego O(1/c) i konserwatywny błąd generatora `<=4A^2/c`.
Dalej opłacono cały rekord próbkowany na siatce oraz zbieżność ścieżek
Skorochoda J1 przez sprzężenie ze zmianą czasu najwyżej 1/c. Uwzględniono
początek, koniec przedziału i więcej niż jeden skok w komórce siatki.

**Nie uzyskano zbieżności pełnych ścieżek w TV.** Podano dokładny wzór
pokazujący jej niezanikający błąd. W dwustanowej kontroli jest to
`1-exp(-1)`, mimo prawidłowej granicy słabej i zbieżności rekordów siatkowych.
Idealne dokładne chwile skoków niosą silniejszą informację niż taka granica;
nie przenoszę jej automatycznie na posterior, filtry lub aparaturę.

## Konkretny warunkowy obiekt aktywności

Dla pierwotnego nieoznaczonego grafu zmian jednej kopii otrzymano
`sqrt(Q_xy Q_yx)=1` na legalnej krawędzi oznaczonych stanów.
Po prawidłowym przejściu do liczności:

    Q(n,n-e_a+e_b)=n_a exp((h_b-h_a)/2).

Zachowano czynniki liczności i stałą przekątną A7 oraz wykazano jawny
limit całkowitej szybkości. Jest to **dostawca aktywności warunkowy względem
nowej hipotezy i jej dopełnienia**, nie niezależne fizyczne źródło HB.
Dotychczasowe prefaktory HB/MET/BAR pozostawiono bez zmian.

## Ważne zastrzeżenie o dopełnieniu

Miara odniesienia K i składnik ukrytych kroków muszą być przekształcane
razem. Jeśli zmienić K w sposób wcześniej równoważny, lecz pozostawić
dopełnienie nieprzekształcone, granica może mieć inne względne aktywności.
To nowy jawny warunek spójności źródła, nie rozwiązanie swobody dopełnienia
przez samą zmianę zapisu. W rzeczywistej badanej rodzinie dopełnienie I
wynika z już zadanej geometrii kroków; nie wyprowadzono jego prawdziwości
fizycznej.

## Weryfikacja i zakres

Przechodzą **54 ograniczone kontrole dokładne**, w tym wymierne przedziały
pierwiastka dla tej samej niejednostajnej pi przy kilku zadanych c,
zewnętrzne przedziały wykładnika, test zmiany czasu, niezgodności dopełnienia
oraz niezależne sumy skoków oznaczonych i liczności. Nie są to 54 dowody
ogólne, symulacja dużego świata ani recenzja innego badacza.

Nośnik obserwowany i pi pozostają stałe; powiększane są matematyczne
jednostajne współrzędne ukryte. Nie jest to rodzimy oddziałujący limit
N->infinity A7 ani theorem z inicjalizacji iid dla niezmienionego HB.
Parametr czasu jest przeskalowanym wewnętrznym indeksem, nie sekundami.

Pełne dowody: [skalowanie](SCALING_AND_GENERATOR_PROOF.md),
[trajektorie i ograniczenie TV](PATH_LIMIT_AND_TV_SCOPE.md),
[aktywność i liczności](ACTIVITY_AND_COUNT_CONSEQUENCE.md).
Strict F01/F03/F06 i fizyczne bramki pozostają otwarte; cel badawczy aktywny.
