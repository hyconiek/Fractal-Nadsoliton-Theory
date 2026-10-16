# KR po audycie: twierdzenie warunkowe zamknięte, uśrednianie pozostaje otwarte

Data:2026-10-16. Szczegóły: [INTAKE_20261016.md](INTAKE_20261016.md).

Kampania wykonała rzeczywisty krok: udowodniła stacjonarną granicę procesu
parzystości dla całej jawnie określonej klasy kinetyk, nie tylko HB i dwóch
wcześniejszych porównań. Warunkiem jest zgodny dokładny dryf, ta sama lokalna
miara Gibbsa, jednokopiowe przejścia, ograniczona aktywność i koncentracja
szumu NC. Słabszy warunek uśredniania po trajektorii IBA też wystarcza.

Przyjąłem ten dowód po związaniu rzeczywistych zależności MT/UF/OBS oraz
sprawdzeniu kontroli reszt, brzegu, maksimów i martyngałów. To twierdzenie
analityczne o modelu, nie formalizacja ani zewnętrzna recenzja.

Sprawdzono46/46 wpisów manifestu. Pięć z sześciu pierwotnych odwołań do
źródeł pasuje; dawny AGENTS pozostaje niedostępny. Odtworzenie producenta
na kopii daje identyczny plik wynikowy. Własny weryfikator wykonuje22
ograniczone kontrole, oddzielając dokładny rachunek od diagnostyki binary64.

Bez NC dowiedziono tylko słabszego wyniku statycznego: średni szum przy
ustalonej współrzędnej wolnej wynosi2 w granicy. Nie znaczy to, że jego
całka wzdłuż każdej losowej historii się uśrednia. Kampania prawidłowo
pozostawiła to pytanie otwarte. Lokalny kolec kowariancji nie jest jego
kontrprzykładem, bo ma zanikającą masę i nadal spełnia NC.

Dodałem [konkretne kryterium dowodowe](REVERSIBLE_POISSON_IBA_CRITERION.md):
oddzielenie niezmiennych składowych i dostatecznie mały koszt równania
Poissona gwarantują IBA. Nie wykazałem jeszcze tego oszacowania dla całej
klasy FIN. Kryterium nie wymaga zgadywania czasu mieszania, ale nie wolno
uznać jego niepowodzenia za dowód zmiany granicy procesu.

Następny łańcuch:
[jawne ograniczenia przepływów i składowych -> oszacowanie Dirichleta/Poissona
-> dowód IBA albo rzeczywisty kontrmodel q12 -> przegląd](../FIN_REVERSIBLE_BRACKET_AVERAGING_NEXT_CHAIN_20261016.md).
To konkretny kolejny problem teoretyczny, bez eksperymentu i nowego prawa
mikroskopowego. Źródło aktywności, naturalne A7, fizyczny zegar, przestrzeń,
zachowana energia, aparatura, QM/grawitacja i ToE pozostają odrębnymi bramkami.
