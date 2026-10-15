# FIN: lokalna równowaga i wykładniczo długie uwięzienie — podsumowanie

## Co nowego

W odróżnieniu od poprzedniego UF (granica trajektorii na każdym *ustalonym* czasie przeskalowanym) celem tej kampanii jest faktyczne okno dla skończonej liczby N kopii: **czy nietknięty proces przez długi czas zachowuje rozkład bliski lokalnej mierze Gibbsa, choć jego pełna równowaga znajduje się gdzie indziej?**

Przedłożony dowód daje dodatnią lokalną barierę `b=5×10⁻¹¹` (celowo konserwatywną, a nie fizyczną wartość), oszacowanie prawdopodobieństwa wyjścia przed czasem T i *wielomianowe*, choć bardzo słabe, oszacowanie szybkości mieszania łańcucha odbitego w tym sektorze.

Wynik syntetyczny w zakresie asymptotycznym:

`sup_{t_loc(N)≤t≤exp(cN)} TV(Law_iid(p_N(t)), pi_D,N) → 0`, `0<c<b`, `N→∞` po wielokrotnościach 12. `t_loc(N)` może mieć absurdalnie wysoki wykładnik wielomianu; nie jest to metoda przewidywania czasu rzeczywistego urządzenia.

## W trzech odrębnych krokach

1. **Uwięzienie:** porównanie z nieznormalizowaną miarą równowagową pozwala szacować przepływ przez wysokoenergetyczną warstwę wyjścia, mimo że start iid-uniform nie jest sam w sobie równowagą.
2. **Jaka lokalna równowaga:** rozkład Gibbsa warunkowany na lokalnym, spójnym obszarze przechodzi w iloczyn rozkładu kwartowego `exp(χ λ z²/2−z⁴/12)` i dziesięciowymiarowego rozkładu Gaussa dla stabilnych fluktuacji zajętości. To NIE jest rozkład pełnej równowagi FIN.
3. **Czy proces dochodzi do tej równowagi:** dla procesu pomocniczego z odbiciem wyprowadzono wielomianowy dolny szacunek luki spektralnej przez ścieżki kanoniczne, konstruowane z rzeczywistych przejść jednej kopii. Rzeczywisty proces i odbity są sprzężone tylko do chwili wyjścia; odrzucenie ścieżek wyjściowych jest ujęte w całkowitym błędzie TV.

## Najważniejsza granica pewności

Całość jest **analitycznym dowodem do niezależnej recenzji**. Niezbędna jest weryfikacja konstrukcji ścieżek, zwłaszcza kroku schodzenia dyskretnego po energii oraz centralnego korytarza `O(log N/N)`. Jeśli ten krok zostanie odrzucony, wolno nadal zachować osobne twierdzenie o uwięzieniu i lokalny rozkład warunkowy, ale nie wolno nazywać okna pełnym dowodem termalizacji.

Opracowanie zakłada istniejącą macierz A7, określony generator leave-one-out, dokładne widmo i matematyczny ciąg `g_N`. Nie wyprowadza żadnego źródła A7 w naturze, fizycznego zegara, cząstki, przestrzeni ani teorii wszystkiego.
