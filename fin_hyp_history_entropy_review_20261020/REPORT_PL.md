# Kontynuacja: zgodność dowodu i odporna niezgodność z HB

Data: 2026-10-20. Drugi przegląd własny z odrębną implementacją kontroli,
nie niezależna recenzja innego badacza. Postulat i oryginalny model FIN
pozostają niezmienione.

## Co sprawdziłem

Przejrzałem całe wyprowadzenie FIN-HYP-HIST-ME-001. W jego zadanym zakresie
utrzymują się istnienie i jednoznaczność prawa, wyprowadzona markowskość,
odwracalność, redukcja do różnic etykiet i niezgodność z HB. Nie znalazłem
potrzeby naprawy matematycznej; nie dopasowywałem nowego postulatu.

Dodałem **drugi dowód istnienia** przez ściśle wypukły, koercywny funkcjonał
dualny. Opłaca on istnienie dodatniego rozwiązania również dla ogólnych
dodatnich danych, nie tylko skonstruowanych przykładów weryfikatora.

## Nowy wynik: małe odstępstwo od maksimum też ma konsekwencje

Dla dowolnej stacjonarnej historii o tej samej mierze stanów i legalnych
zmianach wyprowadziłem

    KL(prawo bloku || prawo maksymalizujące) <= n delta,
    TV <= sqrt(n delta/2),

gdzie delta jest rzeczywistym deficytem entropii na krok, a blok obejmuje
n przejść. Nie zakładam, że porównywana historia jest Markowowska.
Osobno opłaciłem ogon oczekiwania przy przejściu do pierwszego rzeczywistego
wyjścia; nie przenoszę funkcjonału nieograniczonej historii przez TV za darmo.

W istniejącym przykładzie q4 prawo maksymalizujące daje 1/5, a HB 15/71
dla wskazanego pierwszego celu skoku. Każda dopuszczona stacjonarna historia
z tą gałęzią HB musi mieć

    delta >= 21925107363/1739734008083500 > 1/100000

**nata na wewnętrzny krok oznaczonego modelu**. Zatem HB nie odzyskuje się
przez dowolnie małą niedokładność w maksymalizacji. To próg dla konkretnego
wymiernego przykładu, nie liczba A7, wielkość fizyczna lub ostra granica.

## Weryfikacja i punkt zatrzymania

Zachowano wszystkie 18 plików źródła; 17 wpisów manifestu i 11 wejść jest
zgodnych. Odtworzenie 48 kontroli producenta na kopii daje identyczny JSON.
Nowa implementacja przechodzi 68 ograniczonych kontroli, w tym dokładne
przedziały pierwiastka dla nieskonstruowanego p, enumerację krótkich słów
i negatywne testy stacjonarności oraz wyboru miary odniesienia.

Dowody: [przegląd i dual](PROOF_AUDIT_AND_DUAL.md),
[stabilność i przeszkoda HB](ENTROPY_DEFICIT_PATH_THEOREM.md).
Kontrole nie zastępują niezależnej recenzji ani nie dowodzą praw natury.

Ten warunkowy kierunek matematyczny jest dostatecznie zamknięty do decyzji
architektonicznej: hipoteza definiuje inne prawo, a nie źródło dotychczasowego
HB. Strict F01 pozostaje otwarte. Nie rozpocząłem kolejnej odmiany zasady
entropijnej, nie zmieniłem AGENTS, źródeł lub planów, nie wykonałem commita.
