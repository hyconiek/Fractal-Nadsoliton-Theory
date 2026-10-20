# Własna kontrola graniczna AD: rzeczywiste A7 i jawny wynik nierozstrzygnięty

Data: 2026-10-20. Kontrola przedziałowa tego samego autora — NIE niezależna
recenzja, symulacja metastabilności ani nowy dowód całego AD.

Sprawdzono rzeczywiste A7 z odziedziczonych przedziałów spektralnych i
g_c=12/Lambda6, bez podstawiania punktowej wartości zmiennoprzecinkowej.
Cztery z góry zadane populacje obejmują punkt jednostajny, niejednostajny
punkt we wnętrzu, ścianę oraz wierzchołek simpleksu.

Automatyczne różniczkowanie przedziałowe objęło cały zakres
epsilon in [0,1/264] przy każdym z tych ustalonych p. To daje ograniczone
certyfikaty reszty dla nieskończonej sekwencji dopuszczalnych N przy tych
p, nie oszacowanie po całym simpleksie. Bezpośrednio obliczono też wiersze
dla N=264 i 528; liczności są całkowite. Punkty poza centrum testują
globalną część wywodu, nie udają małego-N reżimu MT.

W trzech przypadkach przedziały potwierdzają dodatniość przybliżenia na
całym zakresie epsilon. Przy wierzchołku jest to nierozstrzygnięte przez
ten test: przedział jest zbyt szeroki. Wartości w obu zadanych N są jednak
dodatnie. Nie zmieniono zakresu, punktów ani strojenia, żeby ukryć wynik
niekorzystny. Szeroki przedział reszty przy wierzchołku dodatkowo ostrzega
przed praktyczną interpretacją asymptotycznego twierdzenia.

Warunkowe ograniczenia |F_N-1|<=K_p/N^2 przy ustalonym p są następujące:

| Ustalony punkt | K_p, zaokrąglone na zewnątrz | Dodatniość na całym przedziale epsilon |
|---|---:|---|
| Jednostajny | 3 | potwierdzona |
| Wnętrze niejednostajne | 3 | potwierdzona |
| Ściana | 3 | potwierdzona |
| Wierzchołek | 1450350 | nierozstrzygnięta |

To cztery osobne certyfikaty przy założonych tożsamościach początku
Taylora, nie jedno K działające dla wszystkich p. Wszystkie osiem
bezpośrednich kontroli N=264/528 potwierdza dodatniość w tych punktach.

Dokładny zakres i nierozstrzygnięty punkt są w
[PROOF_OR_OPEN_LEMMA.md](PROOF_OR_OPEN_LEMMA.md) i RESULTS.json. Oryginalny
pakiet AD, AGENTS i wcześniejsze źródła zachowano bez zmian. Samodzielna
kontrola nie jest zamiennikiem
[gotowego zlecenia niezależnego audytu](../fin_global_activity_defect_research_20261020/NEXT_PROOF_AUDIT.md).
F01, fizyczne przesłanki i ToE pozostają otwarte.
Zgodnie z odpowiedzią użytkownika kontynuowano samodzielnie. Dopisano
[analityczną kontrolę argumentu](SAME_AUTHOR_ANALYTICAL_REVIEW.md), w tym
jawny wyraz siatki 4q(theta+h) przy zatrzymanych martyngałach. Nie oznacza
to niezależnej recenzji ani zmiany oryginalnego pakietu.
Końcowy replay obejmuje 39 ograniczonych kontroli, 20 związanych wejść
oraz 331 zachowanych wpisów piętnastu manifestów poprzedników. Nie są to
39 niezależne dowody ogólne.
