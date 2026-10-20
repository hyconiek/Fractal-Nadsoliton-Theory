# FIN — następny krok po badaniu działania KL

**Rezultat:** dokładnie wyznaczono, jakie dane *o historii korelacji* pozwoliłyby jednoznacznie rozpoznać jeden z możliwych generatorów w ścisłej klasie 12-stanowej FIN.

1. Sześć rozróżnialnych modów Fouriera, których stacjonarne korelacje dwuczasowe mierzone są w jednej wspólnej odstępie czasu modelowego, **wystarcza** do odtworzenia sześciu wag przejścia. Macierz odwrotna jest jawna; wyznacznik `-3456 sqrt(3)`.
2. Pięć modów **nie wystarcza**, nawet jeśli dokładnie znamy ich korelacje dla wszystkich chwil. Jawne przekształcenie `w -> w + eta*(-1,1,-1,1,-1,1)` dla dostatecznie małego `eta` pozostawia pięć pierwszych modów niezmienionych, ale zmienia szósty. Wszystkie wagi mogą pozostać dodatnie.
3. Ustalono ograniczenie błędu rekonstrukcji: `||w_hat-w||∞ <= (11/12) eta/[tau(Cmin-eta)]` dla znanej granicy błędów obserwacji `eta`.
4. Sześć wartości korelacji musi ponadto spełniać sześć konkretnych nierówności: `B^-1(-log C)>0`. Przykład `C=(e^-3,e^-1,e^-1,e^-1,e^-1,e^-1)` narusza je mimo dodatniości każdej korelacji. To możliwy **test obalający całą zadeklarowaną klasę modeli**, nie test pochodzenia prawa.
5. **Bez niezależnego zegara** można odzyskać tylko proporcje szybkości, nie fizyczne sekundy.
6. **Bez niezależnego źródła historii** jest to procedura identyfikacji prawa, a nie pochodzenie tego prawa z samego nadsolitonu. F01 nie zostało rozwiązane.

Praca wykorzystuje już istniejącą w teorii matematykę grupy Z12, nie wprowadza nowego fundamentalnego postulatu, bath, zewnętrznej siły ani fizycznej interpretacji stanów jako miejsc w przestrzeni. Nie dotyczy bezpośrednio rank-seven A7, oddzielnego N2 ani pola 1+1.

**Znaczenie dla intuicji:** czasowe korelacje rzeczywiście mogą jednoznacznie kodować wybraną dynamikę w ustalonej klasie. Jednak ich istnienie i wybór nadal wymagają wyjaśnienia. Matematyczne *odczytanie* praw świata z historii nie jest tym samym co ich *wygenerowanie* z informacji.
