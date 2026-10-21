# Strzałka czasu: opłacona entropia odwrócenia rzeczywistych historii FIN

Data: 2026-10-21. Nowy warunkowy dowód analityczny, sprawdzany samodzielnie;
nie niezależna recenzja, formalizacja ani wynik eksperymentalny.

## Wynik

Dla niezmienionego rodzimego HB, zaczynającego od iid i bez odrzucania
niekorzystnych trajektorii, wyprowadziłem zanik względnej entropii między
prawem późnej historii a jej odwróceniem czasowym. Początek historii jest
wybranym, opłaconym czasem wielomianowym tau_N. Wynik jest jednostajny dla
długości zapisu do exp(aN), przy każdym ustalonym 0<a<b z twierdzenia MT.
Obejmuje więc w szczególności okna krytyczne długości sqrt(N)S.

To mocniejsze od dawnego zaniku asymetrii TV. Tamto zastrzeżenie było
poprawne: mały TV sam nie zapewnia małej entropii. Nie usunąłem go przez
założenie stacjonarności rzeczywistego procesu.

## Co było brakującym wyprowadzeniem

Stała mała kula wokół typu jednostajnego leży w prawdziwej składowej D_N.
Daje to wykładniczy ogon wyjątków początkowych iid. Rzeczywiste mieszanie,
oszacowanie wyjść i dokładne podniesienie etykiet płacą porównanie całych
historii. Dodatniość rzeczywistego rozkładu po przygotowaniu ogranicza
logarytm ilorazu wiarygodności przez wielomian. Te dwa oszacowania wolno
już pomnożyć; nie mnożę nieokreślonego o(1) przez niekontrolowany logarytm.

Kluczowe: entropia rzeczywistego prawa względem prawa odbitego jest
nieskończona z powodu różnych nośników. Dowód jej nie używa. Porównuje
rzeczywistą historię z jej własnym odwróceniem; odwracalność generatora
sprowadza iloraz do końców trajektorii.

Rozróżniłem też dwie wielkości: odwrócenie z tym samym przygotowaniem oraz
odwrócony eksperyment zaczynający w rozkładzie końcowym. Dla tego drugiego
modelowa produkcja entropii w natach jest spadkiem D(mu||Pi); jest nie
większa od pierwszej wielkości i również zanika w opłaconym oknie.

W tym konkretnym eksperymencie odwrócenia cała względna entropia pełnej
etykietowanej historii jest już zachowana w dwóch końcowych licznościach.
To dokładny skutek symetrii kopii i brzegowego ilorazu wiarygodności,
nie ogólna wystarczalność dwóch obserwacji dla dowolnych historii czy praw.

## Znaczenie fundamentalne i granice

Te bierne, lokalnie zrównoważone historie nie dostarczają utrzymującej się
strzałki czasu w sensie KL. Dotyczy to także zadanych odczytów parzystości
i całek, z prawidłowym odwróceniem zapisu. Samo gromadzenie wartości
matematycznej ścieżki nie staje się nieodwracalną pamięcią.

Nie wyklucza to nieodwracalnego przygotowania, aktywnej aparatury ani innych
stanów Wszechświata. Nie wyprowadzono fizycznego ciepła, temperatury, zegara
czy kosmologicznej strzałki czasu. Stałe MT, potęga czasu i N0 pozostają
niepraktyczne/nieznane. Wyniku nie przenoszę na każdy początek objęty dawnym
samym o(1) TV. F01 i bramki obserwatora/QM/GR/ToE pozostają otwarte.

[Pełny wywód](PROOF_OR_OPEN_LEMMA.md), [wyniki kontroli](RESULTS.json),
[odtworzenie](REPLAY.md). Znane wzory entropii trajektorii przypisano
literaturze; nowym wkładem jest rzeczywiste opłacenie rodzimego okna FIN.
AGENTS, poprzednie źródła i plany zachowano. Nie dodano nowego prawa.
Przechodzi 45 ograniczonych kontroli; związano 13 wejść i zachowano 30
wpisów dwóch manifestów poprzedników. Nie oznacza to 45 niezależnych
dowodów ani mechanicznego potwierdzenia asymptotycznego twierdzenia MT.
