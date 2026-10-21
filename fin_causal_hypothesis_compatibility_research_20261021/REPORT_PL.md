# Przyczynowość zachowana; dwie hipotezy dają różne następne zdarzenia

Data: 2026-10-21. Nowy dowód warunkowy, samodzielnie sprawdzony;
nie niezależna recenzja ani wyprowadzone prawo natury.

Uwzględniłem wymaganie ciągu przyczynowo-skutkowego. Nie utożsamiam
OPEN_F01 z brakiem dynamiki. Sprawdziłem zamiast tego, czy dwie już
autoryzowane hipotezy rzeczywiście opisują ten sam mechanizm kolejnych
zdarzeń, a nie tylko tę samą statykę.

## Nowy wynik

Dla trzech oddziałujących współrzędnych q12 z zadanym A7 wyprowadziłem
kontrast trzech możliwych następnych zmian liczności. Znika z niego zarówno
wspólny zegar, jak i wolny parametr hipotezy kompozycji kanałów.

Hipoteza kompozycji spełnia tę tożsamość dokładnie. Oryginalna hipoteza
entropii historii narusza ją przy każdym dostatecznie małym dodatnim
sprzężeniu. Znak wynika z dodatniej składowej parzystości A7, nie z dopasowania
liczb. Nie wyznaczyłem użytecznej szerokości tego przedziału i nie przenoszę
wyniku na sprzężenie krytyczne.

To rzeczywista niezgodność następnego zdarzenia, której nie usuwa zmiana
zegara ani wybór parametru drugiej hipotezy. Otrzymałem też ograniczony
kontrast probabilistyczny i dolne oszacowanie odległości TV.

## Znaczenie dla źródła ewolucji

Obie hipotezy nadal definiują przyczynowe prawa modelowe. Ich konflikt
nie oznacza braku przyczynowości. Pokazuje, że nie wolno automatycznie
połączyć dwóch atrakcyjnych zasad i ogłosić, że wspólnie wyprowadzają FIN.
Sama przyczynowość nie wybiera jeszcze jednej z nich; trzeba niezależnie
uzasadnić konkretną własność mechanizmu. Nie dodałem trzeciego postulatu
ani nie wybrałem parametru tak, żeby ukryć konflikt.

## Weryfikacja i przekazanie

[Pełny dowód](PROOF_OR_OPEN_LEMMA.md) i [dokładne współczynniki](EXACT_CERTIFICATE.json).
[Weryfikator](verify_research.py) przechodzi 34 ograniczone kontrole,
w tym 31 dokładnych rachunków oraz 3 diagnostyki binary64. Mały układ
27 stanów sprawdza liniowy rachunek niezależnie od zamkniętego wzoru.
To nie 34 dowody ogólne ani eksperyment fizyczny.

Związano 21 wejść i zachowano 95 wpisów sześciu manifestów poprzedników.
AGENTS oraz wcześniejsze źródła pozostają bez zmian. Nie użyto Firecrawl,
delegacji ani działań zewnętrznych.

[Dalsza decyzja](HANDOFF.md): zachować rozdzielenie obu hipotez i wymagać
niezależnej przesłanki rozstrzygającej mechanizm, zamiast uruchamiać skan
parametrów lub dopisywać nową zasadę naprawczą. Osobny otwarty problem
progowy PDE nie jest przez ten wynik zamknięty. F01, naturalne A7, jednostki,
geometria, aparatura, walidacja, QM/GR i ToE pozostają otwarte.
