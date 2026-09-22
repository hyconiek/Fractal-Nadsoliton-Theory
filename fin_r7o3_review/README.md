# Audyt i integracja R7O3 — Target P

Audyt: 2026-09-21; końcowe przyjęcie: 2026-09-22.
Źródło: `FIN_R7O3_TARGETP_HANDOFF_20260920`.

**Status: PRZYJĘTE — globalne Target P w zadanej rodzinie C4.**
Pełny replay 12 425/12 425 liści oraz końcowe bramki akceptacji przeszły.
Wszystkich 5 432 rodziców zamknięto; nie pozostają nierozstrzygnięte liście.
Wynik zastępuje wcześniejszy status OPEN wyłącznie dla Target P.

Źródłowy katalog pozostawiono niezmieniony. ZIP został porównany z nim podczas
audytu 21 września, ale 22 września nie był już dostępny; nie odtwarzano go.
Końcowa kontrola ponownie sprawdza wszystkie rozpakowane wejścia. Nowe wyniki
kontroli zapisują się tutaj; źródła nie są zastępowane poprawionymi kopiami.

## Dokładny zakres twierdzenia

Dla dostarczonego modelu czterech cech C4 i wspólnych nieujemnych pól
`J3,J4,J5,J6` przyjęto

`lambda2(M4) <= 67/250`, gdzie `M4=Cov_p(C4)`.

Wartości własne są uporządkowane **malejąco**: chodzi o drugą największą.
Równoważnie: najwyżej jedna wartość własna M4 przekracza 67/250.
Stwierdzenie dotyczy zadanej rodziny rozkładów, nie dowolnego rozkładu na
12 etykietach ani pełnej przestrzeni X7.

Konsekwencja dla `H4=I4/g-M4`: najwyżej jeden kierunek o ujemnej krzywiźnie
dla dostarczonego `0<g<=250/67`. Na prawym końcu przedziału samo to ograniczenie
nie wyklucza zerowych wartości własnych. Nie dowodzi dodatniości całego hesjanu.

## Co sprawdzamy i zachowujemy jako dowód

| Element | Zakres |
|---|---:|
| Manifest oraz zgodność katalogu z ZIP-em | 111 plików |
| Oryginalni rodzice resztowi R7N | 5 432 |
| Aktywne certyfikaty R7O3 po świeżym replay | 12 425 / 12 425 PASS |
| Wcześniej przyjęte komórki SAFE R7N | 13 231 |
| Oryginalny podział zwartej dziedziny | 18 663 komórki |
| Podział po rozwinięciu naprawionych rodziców | 25 656 liści |
| Osobna kontrola wymierna | 21 certyfikatów |
| Niezależne kontrole w modelu 12 etykiet | 357 punktów; 2 142 wpisy momentu |
| Testy | 15 nowych + 6 odziedziczonych |

Geometria obejmuje wszystkie drzewa, obie strony każdego podziału, tożsamość
komórek i jednoznaczne przypisanie liści. Nie wystarcza równość sum objętości.
Przyjęte wcześniej 13 231 komórek R7N są identyfikowane z dokładnymi wejściami
ich zakończonego audytu, a nie na podstawie nowych etykiet SAFE w handoffie.

Każdy nowy liść przeliczany jest ze **stałymi zapisanymi wymiernymi B i c**,
bez optymalizatora i bez numerycznego generowania bazy. Wspólne parametry wag
i normalizacja pozostają powiązane w rachunku pochodnych. Pełny replay używa
jawnie zewnętrznie zaokrąglanych przedziałów binary64 z dokładnymi wymiernymi
końcami publicznymi. Końcowy checker dodatkowo odtwarza macierze z nowych
przedziałów momentu i sprawdza ich dodatniość dokładną arytmetyką wymierną.

Osobny proces odtwarza 21 certyfikatów na źródłowym backendzie wymiernym
z siatką `10^-9`, w tym wszystkie 14 źródłowych certyfikatów Gershgorina
oraz liść 10878 o najmniejszym zgłoszonym marginesie. Wszystkie granice
w tej próbie odtwarzają się dokładnie. Nie oznacza to pełnego wymiernego
replay wszystkich liści ani formalizacji w asystencie dowodowym.

Domknięcie dziedziny nieograniczonej wykorzystuje **wcześniej przyjęte FR1
i FR42**. Sprawdzamy ich tożsamość, granice pokrycia i to, że ich próg sigma
jest mniejszy od 67/250; nie liczymy tych dowodów ponownie jako nowych badań.
Odtwarzamy również aktualny provider widma i sprawdzamy normalizację cech.

## Korekty i ograniczenia jakościowe

1. W źródłowym twierdzeniu zapisano rosnący porządek wartości własnych,
   niezgodny z podanym dowodem i konsekwencją. Przyjęty zapis używa drugiej
   **największej** wartości własnej.
2. Zapisany w każdym certyfikacie hash historycznego checkera różni się od
   hasha kodu dostarczonego w paczce. Nie przyjmujemy tego powiązania jako
   potwierdzonego pochodzenia bieżącego kodu. Nowy audyt wiąże świeże wyniki
   z rzeczywiście uruchomionym kodem i jego aktualnym hashem.
3. Źródłowe `verify.py` odczytuje m.in. zapisany raport ukończenia; jego PASS
   nie zastępuje przeliczenia nierówności. Nasza akceptacja wymaga pełnego
   nowego replay, nie odczytu tych flag.
4. Przykład wywołania replay w źródłowym README nie zgadza się z parserem
   argumentów skryptu. Polecenia poniżej są zgodne z nowym checkerem i nie
   usuwają źródłowych checkpointów.
5. Zgłoszone `13/500000000` to margines testu dodatniości, nie globalna
   znormalizowana luka widmowa. Nie wykorzystujemy go do rozszerzenia zakresu g.

Szczegółowy [audyt metody](METHOD_AUDIT.md) oddziela dowód analityczny,
rachunek przedziałowy, kontrolę geometrii, testy punktowe i pochodzenie danych.
Wspólne formuły jetów pozostają współdzielone ze źródłem; niezależny backend
nie oznacza drugiej niezależnej implementacji całego dowodu.

## Czego nie przyjmujemy

Nie domykamy ostrzejszego Target S, pełnego X7, globalnej klasyfikacji minimów
ani ich jednoznaczności. Nie wyprowadzamy źródła gainu, zegara fizycznego,
selektora/QW-2191, dowodów laboratoryjnych, mostu legacy ani transferu ról
fizycznych, SM/GR, `L_total` czy ToE. Istniejące kontrprzykłady X7 pozostają
ważne. Hipotezy z planu kolejnych badań nie stają się wynikami przez to scalenie.

## Artefakty i odtwarzanie

- [Geometria i zależności](registry.json).
- [Świeże przedziały wszystkich przeliczonych liści](leaf_replay.json).
- [Osobny replay wymierny](rational_sample.json).
- [Rekonstrukcja modelu i kontrole analityczne](analytic_checks.json).
- [Pochodzenie i zgodność ZIP-a](provenance.json).
- [Końcowa weryfikacja i przyjęcie](verification.json).
- [Checker](review.py), [końcowe bramki akceptacji](finalize.py),
  [testy](test_review.py).

Z katalogu głównego repozytorium:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o3_review/review.py registry
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o3_review/provenance.py
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o3_review/review.py replay --limit 20 --workers 2
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o3_review/review.py replay --workers 2
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o3_review/review.py replay --rational
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o3_review/analytic_checks.py
PYTHONDONTWRITEBYTECODE=1 python3 fin_r7o3_review/finalize.py --record
```

Pełny replay jest kosztowny i zapisuje atomowy checkpoint co 100 liści.
Zakończony przebieg z dwoma workerami trwał około 18,9 minuty. Dokładna
końcowa kontrola macierzy potwierdziła 12 411 certyfikatów Sylvestera i 14
Gershgorina. Przy braku ZIP-a `provenance.py` zachowuje historyczny raport
porównania i sprawdza dostępne pliki rozpakowane; nie udaje nowego odczytu ZIP-a.
Nie uruchamiaj drugiej jego kopii równocześnie. Nie używaj `python -O`.
Końcowy checker sam nie przelicza ponownie jetów: sprawdza kompletność,
tożsamość danych, dokładne znaki macierzy i testy, korzystając z zakończonego
replay. Nie może awansować niepełnego checkpointu do globalnego PASS.
