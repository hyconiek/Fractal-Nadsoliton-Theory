# Druga hipoteza: składanie wkładów ogranicza kinetykę, lecz jej nie wybiera

Data: 2026-10-21. Wykonane samodzielnie po Twojej zgodzie. To warunkowe
badanie nowej hipotezy, nie przyjęcie prawa FIN ani niezależna recenzja.

## Wynik

Przetestowałem zasadę: po dodaniu dwóch wkładów do potencjału odpowiedzi
szybkości przejścia mają się mnożyć. Przy jawnej symetrii reguły, ciągłości,
niezmienniczości względem stałej potencjału i równowadze Gibbsa otrzymuję
CAŁĄ rodzinę

    R_ab(h)=(1/q)exp[(h_b-h_a)/2+s(h_a+h_b-2 srednia(h))].

Dla q12 pozostaje jeden rzeczywisty parametr s. To rzeczywiste ograniczenie
wcześniej dowolnych aktywności, ale nie ich jednoznaczne źródło. Nie wybrałem
s po obejrzeniu wyniku. Znane szybkości pierwiastkowe są tylko przypadkiem
s=0; nie przedstawiam ich ponownie jako odkrycia.

## Konkretne testy

HB i wcześniejsze przerzedzenie zależne od odległości nie spełniają nowej
zasady w zamrożonym zegarze. Nie znaczy to, że są matematycznie błędne lub
że natura spełnia nową zasadę.

Przy dwóch kopiach s=1/2 może odtworzyć HB po jednej zmianie zegara. Dlatego
sam test N2 byłby mylący. Udowodniłem na rzeczywistym q12/A7, dla N3 i g>0,
że w stanie (0,0,1) dwa konteksty aktualizacji mają ściśle różne normalizacje.
Żadne s nie odtwarza wówczas HB nawet przez wspólne przeskalowanie czasu
zależne od całego aktualnego stanu. Różnicę ujawnia ograniczone zdarzenie:
czy następny skok liczności zmniejsza populację etykiety 1.

W odrębnym dokładnym przykładzie kontrolnym zdarzenie ma prawdopodobieństwa
49/139 i 7/19, z różnicą 42/2641. Te liczby NIE dotyczą A7; dowód dla A7
jest jakościowy i korzysta z jego rzeczywistych własności symetrii/spektrum.
Małe N testuje prawo mikroskopowe, nie domenę metastabilną MT.

Zasada nie gwarantuje jednostkowego zegara prób przy dowolnych polach.
Dla ograniczonego gainu i ustalonego s istnieje natomiast jawny cap CN.
Ponadto liniowa odpowiedź przy stanie jednorodnym i dokładny dryf na linii
parzystości nie zależą od s. Te dane kolektywne nie wybierają aktywności;
nie ogłaszam z tego twierdzenia o granicy trajektorii dla całej rodziny.

## Granica wyprowadzenia i następny krok

Pełny dowód: [PROOF_OR_OPEN_LEMMA.md](PROOF_OR_OPEN_LEMMA.md). Kontrole
rachunkowe i związanie wejść: [RESULTS.json](RESULTS.json). To ograniczone
testy, nie komputerowy dowód wszystkich twierdzeń ani eksperyment.

Wynik tej jednej próby jest terminalny: nowy postulat zawęża kinetykę, ale
nie wybiera s i nie dostarcza rodzimego HB. Nie dodaję teraz kolejnego
aksjomatu, aby naprawić niezgodność. Następny atom źródłowy musiałby wyjaśnić
wpływ wspólnego poziomu potencjału na aktywność, czyli wyznaczyć lub
wykluczyć s z niezależnych przesłanek. Sama znajomość afiniczności,
równowagi albo zgodnego dryfu kolektywnego nie wystarcza.

Naturalne F01 pozostaje otwarte. Nowy postulat nie jest wyprowadzony z
ontologii; A7, Gibbs, nośnik i zegar neutralny są nadal zadane. Poprzednie
hipotezy, źródła, AGENTS i plany zachowano. Nie było commita, delegacji,
wydatków ani przyjęcia fizycznego prawa, przestrzeni, QM, GR lub ToE.
