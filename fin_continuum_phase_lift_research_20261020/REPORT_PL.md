# Jawna realizacja ciągła i odrzucony transfer do skoków FIN

Data: 2026-10-20. Samodzielny dowód analityczny i ograniczone kontrole,
nie formalizacja, niezależna recenzja ani przyjęcie nowego prawa natury.
Model porównawczy: FIN-CONT-PHASE-LIFT-001, nie podstawowy HB ani hipoteza HYP.

## Co skonstruowałem

Istnieje jawny proces odwracalny na tych samych jedenastu zmiennych
granicznych: amplitudzie miękkiej i dziesięciu szybkich. Zachowuje dokładnie
pełną graniczną miarę Gibbsa oraz miękki dryf -z^3/3. Nie dodaje niezależnego
ukrytego bitu — fazę można obliczyć z już istniejącego promienia szybkiego modu.

Realizuje on wcześniejsze profile g_A,g_B jako rzeczywistą wariancję procesu.
Chociaż E[G|Z]=2 w każdym stanie obserwowanej amplitudy, uśrednianie w czasie
zawodzi: oczekiwane maksimum całki z G-2 wynosi dokładnie S M>0.
Zachodzi też różnica dla ograniczonej obserwacji w dwóch chwilach.
[Konstrukcja](PHASE_CHART_AND_GENERATOR.md),
[istnienie i odwracalność](EXISTENCE_REVERSIBILITY_AND_BRACKET.md),
[kontrast zapisu](BOUNDED_RECORD_AND_PATH_CONTRAST.md).

To realizacja MATEMATYCZNEGO kontrmodelu ciągłego. Same statyka i dryf
nie wybierają więc jednego prawa trajektorii w tej klasie. Nie jest to
jeszcze rodzina generatorów na rzeczywistych licznościach FIN.

## Ważny wynik negatywny — zachowany bez poprawiania modelu po teście

Najpierw wyprowadziłem dokładne ograniczenie całkowitego szybkiego nawiasu,
wynikające z każdego pojedynczego skoku liczności. Proste przyspieszenie
części szybkiej konstrukcji przechodzi test rzędu tego śladu, ale to za mało.

Mocniejszy test wariancji konkretnej składowej cosinusowej pokazuje, że
zamrożona konstrukcja narusza lokalny stożek momentów dozwolonych skoków.
Naruszenie występuje na zbiorze o dodatniej masie stacjonarnej, nie tylko
na samej granicy faz. W odpowiedniej normie dwóch charakterystyk pozostaje
dolne ograniczenie błędu 1/32000.

**Ten konkretny transfer charakterystyk do rodzimego FIN został odrzucony.**
Nie zmieniam po tym wyniku współczynników, fazy lub zegara, aby wymusić sukces.
[Ograniczenie skokowe](NATIVE_JUMP_CONE_GATE.md),
[jawny świadek odrzucenia](NATIVE_COSINE_VARIANCE_OBSTRUCTION.md).

Nie wynika stąd niemożliwość każdego generatora zamykającego tę partycję,
ani brak granicy samej amplitudy dla innych szybkich charakterystyk.
Pełne Halla, dokładne szybkości wszystkich liczności, ograniczenie CN
i bezwarunkowe IBA nadal pozostają otwarte. KR pod NC/IBA i nowy wynik MH
zachowują swoje dotychczasowe zakresy.

## Sprawdzenie i dalsza decyzja

Weryfikator przechodzi **47 ograniczonych kontroli**: rachunek skali i
jakobianu, lokalny dual momentów dla wszystkich 132 typów skoków,
rzeczywisty przedział spektralny A7, stałe świadka i osobne małe przykłady
odwracalnego zegara oraz kasowania dryfu. Nie są to 47 dowodów ani
symulacja metastabilności. Związano 13 wejść i zachowano 166 wpisów
dziewięciu manifestów. [Wyniki](RESULTS.json), [replay](REPLAY.md).

Przez browser-harness przeczytałem fragmenty pierwotnej pracy Li–Ying
o symetrycznych dyfuzjach. Jej dowolne przybliżenia łańcuchowe nie opłacają
ograniczeń konkretnego grafu FIN. [Zakres lektury](LITERATURE_NOTE.md).

Ten etap kończy się obiektem ciągłym ORAZ konkretnym wynikiem przeciwnym
próbie jego natywnego transferu. Nie uruchamiam kolejnej serii wariantów.
F01, naturalne źródło kinetyki, jednostki, geometria, fizyczna walidacja
i QM/GR/ToE pozostają otwarte. Źródła, plany i AGENTS zachowano.
