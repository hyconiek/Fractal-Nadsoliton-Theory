# Audyt i scalenie badań rank-7 — 2026-09-19

Ten katalog jest teraz wspólnym punktem wejścia do zweryfikowanych wyników
pakietów z 14–16 września. Oryginalne katalogi pozostają niezmienionymi
archiwami; ich nagłówki `CERTIFIED`, `DONE` i `PASS` nie zastępują tego audytu.

## Wynik

- Potwierdzono brzegowe ograniczenie Isinga, ograniczenie części
  wewnątrz klas parzystości oraz lokalne wyniki R7P w opisanych zakresach.
- Potwierdzono FR1 i silniejszy zakres FR42: `exp(-2 J6)<=1/1000000`.
  Dla FR42 ponownie sprawdzono wszystkie 637 liści i lokalne przesłanki.
- Z 106 prostokątów FR223 **99 przechodzi bez podziału**. Pozostałe siedem
  nie przechodziło dostarczonego testu na całym prostokącie. Naprawiono je
  przez podział na łącznie **18 certyfikowanych podprostokątów**. Dopiero te
  dodatkowe świadectwa uzasadniają przyjęcie całej niepowiększonej unii.
- Ponownie certyfikowano **co najmniej 60 lokalnych pierwiastków kwartowych
  i 60 pełnych**, dla dokładnie określonych dziesiętnych amplitud. Zastąpiono
  zaokrąglone stałe prawdziwymi przedziałami π i pierwiastków oraz sprawdzono
  bezwładność hesjanu metodą przedziałową. Nie jest to wyczerpanie wszystkich
  pierwiastków ani dowód globalnej równoważności obu krajobrazów.
- **119 testów w 29 plikach przeszło** na świeżej kopii dla każdego pliku.
  Dodatkowo wykonano opisane niżej kontrole naukowe, które nie są tylko
  odczytem zapisanych flag `PASS`.
- Przeszło także **7 nowych testów integracyjnych**, w tym odrzucenie
  usuniętego liścia pokrycia, fałszywego znaku, zmienionej geometrii i błędnego
  punktu startowego certyfikacji fazowej.

Globalne ograniczenie czterowymiarowe, pełne wyczerpanie punktów stacjonarnych
7D i jednoznaczność globalnego minimum nadal pozostają otwarte.

## Czy nowsze katalogi zawierają starsze?

Porównanie dotyczy plików roboczych, z pominięciem cache Pythona/pytest.

| Porównanie | Starszych plików | Identyczne | Zmienione | Brakujące | Dodane |
|---|---:|---:|---:|---:|---:|
| starszy bundle → `fin_rank7_followup` | 150 | 141 | 9 | 0 | 140 |
| `fin_rank7_followup` → kontynuacja 20260915 | 290 | 285 | 4 | 1 | 29 |
| kontynuacja 20260915 → FR223 | 318 | 314 | 4 | 0 | 52 |

Jedyny brak między 14 a 15 września to plik pomocniczy
`tests/test_G_boundary_cover.py.tmp`, nie brak wyników naukowych.
Nowszy FR223 zawiera wszystkie względne ścieżki pakietu z 15 września,
ale cztery dokumenty/ledgery zostały zaktualizowane. Nie jest więc identyczną
kopią bajtową całego starszego katalogu. Zewnętrzny historyczny handoff
`FIN_rank7_followup_handoff_current.md` i manifest opakowania starszego bundle
nadal należy zachować jako historię pochodzenia.

Manifesty dostarczonych kontynuacji przeszły kontrolę. Stare manifesty wewnątrz
odziedziczonych pakietów opisują wcześniejsze checkpointy i nie są aktualnymi
manifestami końcowymi. [Pełne porównanie](continuation_inventory.json).

## Co poprawiono metodologicznie

1. **Testy modyfikujące pakiet.** Wcześniejszy sekwencyjny replay w tym katalogu
   dał 93 sukcesy i jedną porażkę: test zapisywał na nowo certyfikat R7P-068
   bez wymaganych metadanych, co psuło późniejszy test weryfikatora. To błąd
   izolacji testów, nie obalenie matematyki. Nowy replay używa świeżych kopii
   w `/tmp`. Oryginały i wcześniejszy zapis porażki zostały zachowane.
2. **Brak pytest w lokalnym interpreterze.** Użyto jawnego runnera `unittest`
   z obsługą bezargumentowych funkcji testowych; nie udajemy uruchomienia
   pytest. Wszystkie 119 dostarczonych testów zostało wykonanych, nie pominiętych.
3. **Weryfikator metadanych nie jest weryfikatorem dowodu.** Dostarczony
   `verify.py` kontroluje głównie hashe i schematy. Część testów czyta jedynie
   zapisane wyniki. Dlatego osobno sprawdzono obliczenia, dowody i pokrycia.
4. **FR223: granice pojedynczych prostokątów.** FR32, FR48, FR52, FR54, FR56,
   FR58 i FR60 wymagają nowych świadectw podziału. Nie uznano ich za poprawne
   tylko dlatego, że występowały w masce wyszukiwania. Pozostałe 99 obszarów
   nie pokrywa w całości tych siedmiu; sprawdzono to osobno.
5. **Bufor 1.02.** Jest wyłącznie narzędziem nawigacji optymalizatora. Nie
   powiększa udowodnionych obszarów. Brak dodatniej luki w czterech próbach
   FR223 nie dowodzi braku kontrprzykładu w pozostałej dziedzinie.
6. **Luka macierzy Schura.** Odrzucony współczynnik `2/25` pozostaje odrzucony.
   Bezwładność przy zadanym progu nie przenosi liczbowych odległości wartości
   własnych między macierzą Schura a M4.
7. **Stałe fazowe i bezwładność.** Pierwotne checkery używały m.in.
   `math.pi`, `math.sqrt` i numerycznych wartości własnych bez osobnego
   certyfikatu błędu. Nowy checker używa przedziałowych stałych, ścisłej
   inkluzji Krawczyka, warunku kontrakcji i przedziałowego LDL po odwracalnej
   zmianie bazy. Dotyczy to ustalonego fixture, nie amplitud nieznanego dokładnego
   punktu koegzystencji.
8. **Średnica C4.** Odległości cech C4 nie muszą zależeć tylko od odległości
   cyklicznej etykiet. Poprawiono uzasadnienie pomocnicze FR42: sprawdzono
   wszystkie 66 par, zachowując wystarczające ograniczenie `D²<5`.
9. **Historyczne `delta_L`.** Różnica `5e-16` wynikała z utraty zależności przy
   odejmowaniu oddzielnie zaokrąglonych przedziałów Laplasjanu. Oryginalny
   provider widma W odtwarza dokładnie historyczne ograniczenie. Nie znaleziono
   tu błędu wcześniejszego twierdzenia o discord.
10. **Niepełne świadectwa.** FR2/FR3 i część historycznych kroków FR4–FR7 oraz
    FR20/FR21 nie są w tym audycie awansowane do odtwarzalnych certyfikatów.
    Brakuje ich pełnych generatorów/pokryć. Nie użyto ich jako przesłanek
    nowych zaakceptowanych wyników. Brakuje też wskazanego w starszym opisie
    pliku `R7P-077_local_jacobians.json`; precyzyjnego certyfikatu alignment
    R7P-079 nie przyjęto wyłącznie na podstawie jego zapisanego podsumowania.

## Przyjęte wyniki i kod

[ACCEPTED_RESULTS.md](ACCEPTED_RESULTS.md) podaje dokładne zakresy naukowe
oraz zależności. [verify_continuation.py](verify_continuation.py) sprawdza
geometrię, pokrycie FR42 i prostokąty. [scientific_rechecks.py](scientific_rechecks.py)
zawiera poprawione certyfikacje fazowe, świadka energetycznego i punktu
niestabilności kątowej. [review.py](review.py) wykonuje porównania i izolowane
testy; [run_test_file.py](run_test_file.py) jest jawnym runnerem kompatybilności.

Nowe dane kontroli:

- [119 testów](replay_latest_fresh.json).
- [99 bezpośrednich certyfikatów / 7 początkowych odrzuceń](FR223_union_replay.json).
- [18 podobszarów naprawczych](FR223_subdivision_repairs.json).
- [Pokrycie FR42](FR42_replay.json) i [niezależna geometria FR1](FR1_independent_geometry.json).
- [120 poprawionych certyfikatów lokalnych](phase_recertification.json).
- [Energia, kąt i odtworzone delta_L](global_rechecks.json).

Do obliczeń prostokątowych dodano zewnętrzne zaokrąglanie wymierne do siatki
`10^-60` po każdej operacji. Zaokrąglanie odbywa się zawsze na zewnątrz,
wyłącznie poszerzając przedziały. Przyspiesza to kontrolę bez zastępowania
rachunku przedziałowego obliczeniem zmiennoprzecinkowym. Zachowano też
częściowy replay bez tego uproszczenia, który potwierdził m.in. odrzucenie FR32.

## Odtworzenie

Uruchamiaj z głównego katalogu repozytorium. Wszystkie polecenia zapisują nowe
ledgery wyłącznie tutaj; źródła badań pozostają niezmienione.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/review.py --latest-fresh
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/verify_continuation.py inventory
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/verify_continuation.py geometry
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/verify_continuation.py tails
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/verify_continuation.py boxes
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/verify_continuation.py repair_boxes
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/scientific_rechecks.py phases
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/scientific_rechecks.py globals_recheck
PYTHONDONTWRITEBYTECODE=1 python3 fin_rank7_intake_review/finalize.py
```

Ostatnie polecenie sprawdza integralność i strukturę zapisanych świadectw oraz
uruchamia siedem testów integracyjnych; nie zastępuje pełnego naukowego replay.
[Końcowy rejestr](consolidated_verification.json) zawiera hashe i wyniki.

Hashe i istniejące logi nie zastępują ponownego obliczenia dowodu. Nie wykonano
pełnego globalnego pokrycia pozostałej dziedziny ani nowej kampanii FR224.
Nie wykonano commitów, nie usunięto starszych pakietów i nie wygenerowano PDF.
