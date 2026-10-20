# FIN — plan wykonania i decyzji

Data: 2026-10-19. **Plan niewykonany.** Architekt przygotował zadania;
nie uruchomił kampanii ID, RP004 ani MH.

## 1. Domyślny start

Uruchomić wyłącznie [ID-001](research_tasks/01_ID001_ENTROPY_SOURCE_BOUND.md).
Sol nie rekonstruuje całej historii: korzysta z mapy i wskazanych wejść,
wykonując tylko ukierunkowaną kontrolę świeżości danego celu.

Gotowe polecenie:

```text
Pracuj w /media/magazyn/home/teresa/Dokumenty/HYC/TOE/edison.
Przeczytaj aktualne AGENTS.md, research_tasks/00_EXECUTION_RULES.md
i research_tasks/01_ID001_ENTROPY_SOURCE_BOUND.md.
Wykonaj tylko ID-001 w jego nowym katalogu wynikowym.
Cele dowodowe nie są założeniami ani przyjętymi wynikami.
Wykorzystaj wskazane poprzedniki i nie powtarzaj audytu całego repo.
Zwróć dowód, poprawny kontrprzykład lub dokładny open lemma,
niezależną kontrolę kluczowego kroku i kompletny handoff.
Nie zmieniaj AGENTS ani praw FIN. Po wyniku zatrzymaj się.
```

## 2. Decyzja po ID-001

| Wynik | Następna czynność |
|---|---|
| Budżet i native corollary opłacone | Astra ocenia wykluczoną klasę źródeł i rzeczywiście dostępne alternatywy. Bez nowego obiektu źródło pozostaje OPEN_INPUT. |
| Tylko finite-delta inequality | Zachować ją; nie głosić all-resolution no-go bez opłaconych kwantyfikatorów. |
| Nieopłacone finite-H w aplikacji | Ograniczyć zastosowanie; nie dopisywać finite-H do ontologii nadsolitonu. |
| Kontrprzykład | Wskazać upadłą implikację, sprawdzić miarę, joining, zegar i definicję entropii; wrócić do Astry. |
| Późniejszy równoważny wynik | Reuse i DUPLICATE_TARGET, bez drugiej kampanii. |

Znane ogólne infinite-entropy/no-exact-factor nie jest nowym wynikiem FIN.
Wkładem zadania ma być konkretny kontrakt błędu, rozdzielczości i źródła.
Sukces tego ograniczenia nie wyprowadza sam nowego prawa zmiany.

## 3. Osobno wybrana kontynuacja modelu polowego

[RP004-A](research_tasks/02_RP004_A_RETARDED_REFERENCE.md) i
[RP004-B](research_tasks/03_RP004_B_UNIFORM_REMAINDER.md) mogą pracować
równolegle po wybraniu tej gałęzi. Czytają ten sam
`research_tasks/RP004_SHARED_CONTRACT.json` i mają rozłączne katalogi zapisu.
Nie zmieniają kontraktu ani wyników drugiego zadania. Nie więcej niż dwa
aktywne pakiety; równoległość dowodowa nie jest zgodą na ciężkie symulacje.

[RP004-C](research_tasks/04_RP004_C_SIGNED_FLUX.md) wymaga obu naukowo
opłaconych oszacowań. Nie wystarczy obecność plików lub status PASS.

| A/B | Reguła dla C |
|---|---|
| Retarded lower bound i trace error w tym samym oknie | Przeliczyć kompletny oryginalny podpisany bilans. |
| A tylko fixed-R | Brak opłaconego R=T/6; zwrócić konkretny brak. |
| B tylko fixed-T lub H1 x L2 | Nie zakładać jednostajności lub śladów; zwrócić brakujący bound. |
| Niezgodne dane, normalizacje, okna lub zakresy parametrów | Powrót do Astry przed łączeniem. |
| Kontrprzykład do zamrożonego celu | Zamknąć ten cel w jego zakresie, zachować NS i Gamma. |

C obejmuje również terminalny handoff istniejącego RP004-D. Nie tworzymy
osobnej kampanii wyłącznie raportowej. Po wyniku C wrócić do architekta,
bez automatycznego badania absorpcji, nowych kerneli lub interpretacji fotonowej.

## 4. Opcjonalny Hall

[MH-001](research_tasks/90_MH001_MESOSCOPIC_HALL.md) jest matematycznie
niezależny od ID/RP, ale ma niższy priorytet fundamentalny. Jedna kula,
L=1/10, SHELL i ADJACENT. Sukces dotyczy tej rodziny cięć, nie pełnego Hall,
globalnego gluing, cap5N ani IBA counterfamily. Porażka wystarczającej metody
nie dowodzi uniwersalnego IBA. Przed zmianą rodziny wrócić do Astry.

## 5. Poza kolejką

- AS/KB: bez rzeczywistego niezależnego obiektu nie uruchamiać pustej bramki.
- Źródło A7, nieliniowości, jednostek i selektora: nowe wejście lub nowy
  blocker-cut są konieczne; stare inwentarze nie stanowią nowej pracy.
- B4: odroczyć do celu wymagającego globalnego wykładnika i pełnego replay.
  Task335 nie rusza przed zamrożeniem właściwego predyktora.
- Empiryka: rzeczywiste dane, kalibracja, zasoby i osobna autoryzacja;
  nie generować symulacji nazwanej walidacją.
- Nowy postulat: osobny model i decyzja architektoniczna przed przyjęciem.

To są bramki, nie fikcyjnie gotowe zadania z brakującymi wejściami.

## 6. Sol i Astra

Sol wykonuje zamrożone lematy, kontrolne rachunki i kontrprzykłady, wiąże
źródła, naprawia lokalny błąd w ramach metody i oddziela dowód od numeryki.
Nie potrzebuje decyzji architekta przy każdym rachunku.

Astra wraca przy nowym źródle dynamiki, sprzeczności z przyjętym twierdzeniem,
zmianie klasy dopuszczalnych mechanizmów, nośnika, obserwacji lub zegara,
niezgodności A/B, pełnym theorem/counterfamily IBA, albo propozycji fizycznego
interfejsu. Każdy pakiet dopuszcza jeden wskazany kierunek dowodu i najwyżej
jedną naprawę wynikającą z konkretnego błędu; potem dokładny open lemma i stop.

## 7. Weryfikacja

Sprawdzać zakres, źródła, wszystkie krytyczne zależności, znaki i normalizacje,
niezależną kontrolę wrażliwego kroku, zachowanie negatywnych wyników i brak
nieopłaconej promocji fizycznej. Nie powtarzać kosztownych niezależnych
kampanii tylko dlatego, że istnieje nowy raport.

Zachować źródła, AGENTS i historyczne hashe. Nowy live context wiązać osobno.
Nie wznawiać zatrzymanej ciężkiej algebry. Bez commit/push i działań zewnętrznych.

Maszynowy plan: [research_tasks/TASKS.json](research_tasks/TASKS.json).

```bash
python3 fin_strategy_audit_20261019/verify_strategy.py
```

Ten checker sprawdza architekturę, źródła i zależności. Nie dowodzi żadnego
z przyszłych twierdzeń. Wszystkie cele ID/RP/MH pozostają niewykonane.
