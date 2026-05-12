# P1313 — O3 Nonuniqueness Exclusion (strict-only) — execution plan, theorem skeleton, and audit protocol (PL)

Status: `O3_IN_PROGRESS_NOT_CLOSED`
As of: `2026-05-12`
Depends on: `P1312`

## Cel dokumentu
Ten pakiet realizuje **następny uczciwy krok z P1312**: uruchomienie `O3 Nonuniqueness exclusion` jako bramki krytycznej przed jakimkolwiek claimem `QW-2191 = CLOSED` w lane strict-only.

Dokument:
1. nie ogłasza domknięcia `QW-2191`,
2. nie miesza strict-only z NB-only,
3. nie używa mostu legacy->strict,
4. definiuje formalny plan dowodu + protokół kontrprzykładów + warunek przejścia `nonuniqueness_residual: OPEN -> CLOSED`.

## Guardrail compliance (jawnie)
- `K_strict_gate` pozostaje operacyjnym obiektem strict lane; brak transferu ról legacy bez osobnego twierdzenia.
- `QW-2191` traktowane jako realna przeszkoda unikalności, nie jako "prawie zamknięte".
- Brak implikacji: `NB closed => strict closed`.

## O3 — sformułowanie formalne (wersja robocza)

### Definicje robocze
Niech `C_strict` oznacza klasę dopuszczalnych strict-core kandydatów selektora (zgodnych z aktualnymi guardrails i bez axiom-tagged domknięć non-strict).

Niech `Dir(c)` będzie kierunkiem indukowanym przez kandydata `c in C_strict` po translacji do wspólnej reprezentacji porównawczej.

Niech `~` oznacza relację równoważności operacyjnej kandydatów (ten sam kierunek + ten sam wynik na zbiorze testów zgodności strict).

### Twierdzenie O3-EXCLUSION (target)
Dla każdego `c1, c2 in C_strict` zachodzi dokładnie jedno:
1. `c1 ~ c2` (równoważność), albo
2. co najmniej jeden z kandydatów jest niedopuszczalny (narusza constraints strict-core),

a ponadto nie istnieje para dopuszczalnych kandydatów strict-core taka, że `Dir(c1) != Dir(c2)`.

### Konsekwencja
Jeżeli twierdzenie O3-EXCLUSION zostaje udowodnione i replay/adversarial audit przejdą, wtedy:
- `nonuniqueness_residual = CLOSED`.

Jeżeli którykolwiek warunek zawiedzie:
- `nonuniqueness_residual = OPEN`.

## Plan dowodu (profesorski, minimalny i audytowalny)

### Krok O3.1 — Enumeracja klasy kandydatów
- Jawnie wypisać skończony katalog klas `C_strict = {C1, ..., Cn}`.
- Dla każdej klasy podać:
  - interfejs wejście/wyjście,
  - ograniczenia dopuszczalności,
  - mapę translacji do wspólnej reprezentacji `Dir(.)`.

**Warunek PASS O3.1:** brak "ukrytych" klas i brak otwartych aliasów parametrycznych.

### Krok O3.2 — Lemat normalizacji reprezentacji
- Udowodnić, że translacja `T: C_strict -> R_common` jest jednoznaczna w granicach tolerancji strict.
- Wykluczyć gauge-like redundancje, które mogłyby tworzyć sztuczną wieloznaczność.

**Warunek PASS O3.2:** każdy kandydat ma pojedynczy obraz w `R_common`.

### Krok O3.3 — Lemat rozdzielający kierunki
- Dla każdej pary klas `(Ci, Cj)` przeprowadzić test:
  - albo wykazać `Dir(Ci) = Dir(Cj)` po normalizacji,
  - albo wskazać formalny powód niedopuszczalności jednej klasy.

**Warunek PASS O3.3:** brak pary dopuszczalnych klas z rozbieżnym kierunkiem.

### Krok O3.4 — Counterexample sweep
- Uruchomić checklistę kontrprzykładów:
  1. perturbacja fazowa,
  2. perturbacja damping/exponent w granicach strict tolerance,
  3. edge-case na granicy admissibility,
  4. próba wymuszenia alternatywnego kierunku przez ukrytą reparametryzację.

**Warunek PASS O3.4:** `counterexample_count = 0` dla dopuszczalnych klas.

### Krok O3.5 — Replay + adversarial audit
- Niezależny replay przebiegu O3.1–O3.4 w innej kolejności testów.
- Adversarial próbuje złamać unikalność bez naruszenia strict constraints.

**Warunek PASS O3.5:** replay = PASS, divergence = 0, adversarial_break = FAIL.

## Macierz decyzji O3
- `O3 = PASS` tylko gdy `O3.1..O3.5 = PASS`.
- W przeciwnym razie `O3 = FAIL` i `nonuniqueness_residual` zostaje `OPEN`.

## Status na dziś (uczciwy)
- O3.1: `NOT_EXECUTED`
- O3.2: `NOT_EXECUTED`
- O3.3: `NOT_EXECUTED`
- O3.4: `NOT_EXECUTED`
- O3.5: `NOT_EXECUTED`
- Globalnie: `O3_IN_PROGRESS_NOT_CLOSED`

## Czego ten dokument nie twierdzi
- Nie twierdzi, że `QW-2191` jest już domknięte.
- Nie twierdzi, że istnieje strict-core selector closure już teraz.
- Nie twierdzi globalnego ToE closure.

## Rekomendowany następny uczciwy krok
Wykonać **O3.1 Enumerację klasy kandydatów** i zafiksować skończony katalog `C_strict` z jawnymi kryteriami dopuszczalności. Bez tego dalsze testy O3.2–O3.5 są metodologicznie niedomknięte.

## Dla laika
To jest etap "spisu wszystkich legalnych kompasów" zanim sprawdzimy, czy każdy pokazuje ten sam kierunek. Jeśli nie mamy pełnej listy legalnych kompasów, nie da się uczciwie udowodnić, że kierunek jest unikalny.
