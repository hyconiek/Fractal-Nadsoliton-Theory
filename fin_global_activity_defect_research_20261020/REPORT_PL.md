# Granica krytyczna oryginalnego źródła: opłacono brakujące oszacowania

Data: 2026-10-20. Nowy dowód warunkowy, samodzielnie audytowany; nie
niezależna recenzja, formalizacja ani przyjęcie prawa natury.

## Wynik

Wyprowadziłem rzeczywisty krytyczny dryf i wariancję **oryginalnego źródła
hipotezy na N oddziałujących współrzędnych**. Daje ono granicę trajektorii

    dZ=12(chi lambda Z-Z^3/3)ds+sqrt(24)dB,
    Z0~rho_lambda,

w ustalonej wcześniej normalizacji wewnętrznych kroków. Nie dodano
jednostajnych współrzędnych, kąpieli ani nowego prawa szumu. Wynik obejmuje
lokalny start stacjonarny oraz opłacone przejście do późnego startu iid.

## Co było rzeczywiście nową trudnością

Samo przybliżenie normalizacji nie wystarczało. Wyprowadziłem poprawkę
pierwszego rzędu, która daje błąd marginesów O(N^-2), oraz lokalną
tożsamość porównującą dodatnie przepływy z dokładnym globalnym rozwiązaniem.
Zachowanie wyrazu brzegowego usuwa problem wykładniczo małej masy lokalnego
sektora w pełnym Gibbsie. Ważone oszacowanie pozwala kontrolować dryf,
nie tylko średnią macierz lub nieprzeskalowany błąd prawdopodobieństw.

Nie pominięto też dyskretnej poprawki „wariancja = drugi moment minus
kwadrat średniej”. Kontrola martyngałów, maksimum, wspólnej filtracji i
rzeczywistych początkowych wyjątków jest częścią twierdzenia.

## Czego to nie zamyka

Jest to konsekwencja upoważnionej hipotezy i zadanych A7/Gibbsa/nośnika,
nie ich naturalne źródło. F01/F03/F06, fizyczne jednostki, czas, geometria,
aparatura, walidacja, QM/GR i ToE pozostają otwarte. Nie udowodniono pełnej
zgodności historii z HB/SQRT ani stałej punktowej aktywności na każdej
krawędzi. Starych negatywnych wyników nie usunięto.

[Pełny dowód](PROOF_OR_OPEN_LEMMA.md) wymaga teraz
[niezależnego audytu](NEXT_PROOF_AUDIT.md); nie proponuję kolejnej serii
parametrów lub cyfr. Źródła i AGENTS zachowano. Kontrole w RESULTS.json
sprawdzają ograniczony rachunek, nie zastępują dowodu asymptotycznego.
Przechodzi 39 ograniczonych kontroli dokładnych. Związano 15 wejść oraz
zachowano 311 wpisów czternastu manifestów poprzedników. Początkowy problem
zapisu bardzo długiego ułamka i jego naprawę zachowano w notatce weryfikacji;
nie zmieniał on żadnej nierówności naukowej.
