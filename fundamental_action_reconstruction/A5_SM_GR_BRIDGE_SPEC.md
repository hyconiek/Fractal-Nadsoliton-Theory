# A5 Spinor Route Split

Status: `A5_EXECUTED_SPINOR_ROUTE_SPLIT_WITH_REPO_PRIOR_ART`
As of: `2026-03-06`

## Uwaga o nazwie pliku

Nazwa pliku `A5_SM_GR_BRIDGE_SPEC.md` zostaje utrzymana dla kompatybilnosci repo.
Rzeczywisty zakres etapu `A5` jest teraz wezszy i poprawniejszy:
- najpierw rozdzial drogi spinorowej,
- dopiero potem wejscie do `A6` i dalszych mostow gauge / gravity.

## Cel

Po `A1-A4` rozdzielic dwie dopuszczalne drogi dla sektora fermionowego:
1. `spinor-emergent route`
   - zgodna z ontologia jednego nadsolitonu,
   - wymagajaca 3D topologii, Hopf/Skyrmion structure i analizy operatorowej.
2. `minimal spin-bundle extension`
   - mniej ambitna ontologicznie,
   - ale szybsza jako kontrolowany most efektywny do Dirac/gamma sector.

To jest etap klasyfikacji i decyzji konstrukcyjnej.
Nie jest to jeszcze theorem-level derivation spinorow.

## Wczesniejsze badania w repo relewantne dla A5

Nie wszystkie poprzednie badania sa metodologicznie rowne.
Dlatego `A5` nie traktuje calego legacy korpusu jako dowodu.
Wprowadzone zostaje twarde rozroznienie:
- `strict-admissible internal references`
- `legacy / exploratory / methodology-risk`

Tylko pierwsza grupa moze byc uzywana jako wewnetrzny material konstrukcyjny.
Druga grupa moze byc uzywana najwyzej jako:
- heurystyka,
- generator hipotez,
- negatywna kontrola,
- sygnal, czego nie wolno jeszcze sprzedawac jako domknietego.

### Strict-admissible internal references

1. `QW-2121`
   - istnieje formalna kompletna specyfikacja rozszerzenia spinor+gauge,
   - ale status pozostaje `formal_extension_not_strict_derived`.
2. `QW-2126`
   - istnieje deterministyczny numeric bridge dla `g`, `g'`, `g3` i Yukaw z pakietu strict,
   - bez scan/no-fit.
3. `QW-2127`
   - istnieje action-level nonabelian spinor+gauge bridge,
   - ale `representation_assignment_unique_from_kernel=False`,
   - `anomaly_cancellation_proof_from_kernel=False`.

### Warstwy de-anchored i kernel-mode

4. `QW-2189`
   - de-anchored spinor+gauge consistency:
   - exact anomaly/charge closure,
   - ale full representation emergence pozostaje otwarte.
5. `QW-2190`
   - istnieje kernel-mode representation scaffold z embedded Lie-closure `SU(3)xSU(2)xU(1)`,
   - ale fizyczna unikalnosc mapowania nie jest zamknieta.
6. `QW-2191`
   - istnieje jawny obstruction theorem:
   - degenerate eigenspaces daja ciagla rodzine `O(2)` i blokuja pelna axiom-free unikalnosc mapowania.

### Legacy / exploratory / methodology-risk

7. `QW-1200`
   - 3D skyrmion / Hopf / Finkelstein-Rubinstein / Jackiw-Rebbi route daje wiarygodny mechanizm fermionowego emergence,
   - ale pozostaje poza theorem-level closure,
   - i nie jest tutaj traktowany jako dowod, tylko jako kierunek heurystyczny.
8. `20 HYDRODYNAMIC PROPERTIES OF INFORMATION SUPERSOLITON...`
   - jawnie wykazano, ze prosty 2D vortex route nie daje spinu `1/2`,
   - ale i ten wynik jest traktowany jako legacy signal / negatywna kontrola, nie jako finalny proof artifact.

## Route 1: Spinor-emergent

### Idea

Spinor nie jest doklejany recznie.
Ma wyjsc z topologii i spektrum jednego nadsolitonu po przejsciu do poprawnej 3D struktury.

### Minimalny pakiet potrzebny

1. 3D topologiczny obiekt typu Skyrmion/Hopf z ladunkiem jednostkowym,
2. mechanizm Finkelstein-Rubinstein dla znaku wymiany,
3. Jackiw-Rebbi-like zero modes lub rownowazna analiza operatora liniaryzacji,
4. wyprowadzenie Dirac-like effective operatora z fluktuacji wokol tla nadsolitonowego,
5. kompatybilnosc z kernel-mode scaffold `QW-2190`.

### Co juz wspiera ten route

- `QW-2189/2190` daja warstwe consistency + partial representation scaffold,
- ontologicznie ten route najlepiej pasuje do `single nadsoliton`.

### Co wolno powiedziec uczciwie

- 3D topological route jest zgodny z ontologia teorii,
- ale jego bezposrednie wsparcie w repo pochodzi czesciowo z badan legacy o niepewnej metodologii,
- wiec route jest traktowany jako `primary hypothesis branch`, a nie jako lokalnie potwierdzony derivation package.

### Co pozostaje otwarte

- brak theorem-level derivation operatora Diraca z `A1-A4`,
- brak jawnego wyprowadzenia macierzy gamma z liniaryzacji nadsolitonu,
- brak pelnej unikalnosci reprezentacji przez obstruction `QW-2191`,
- brak zamkniecia anomalii bez residual external template assumptions na poziomie czysto kernelowym.

### Ocena

- zgodnosc z ontologia: `wysoka`
- zgodnosc z admissible strict references: `czesciowa`
- gotowosc wykonawcza: `niska do sredniej`
- rygor recenzencki dzis: `partial only`

## Route 2: Minimal spin-bundle extension

### Idea

Spinor jest dolaczany jako kontrolowane minimalne rozszerzenie geometryczne nad efektywna warstwa metryczna, bez twierdzenia, ze sam juz emerguje z nadsolitonu.

### Minimalny pakiet potrzebny

1. efektywna metryka / tetrada / spin connection,
2. standardowy blok Dirac:
   - `bar(psi) (i gamma^mu D_mu - m) psi`,
3. zgodnosc z action-level bridge `QW-2121/2126/2127`,
4. zachowanie kompatybilnosci z anomaly/charge closure `QW-2189`.

### Co juz wspiera ten route

- `QW-2121`: kompletna formalna specyfikacja rozszerzenia,
- `QW-2126`: numeric bridge dla couplingow i Yukaw,
- `QW-2127`: action-level nonabelian bridge,
- `QW-2189`: de-anchored consistency layer.

### Co pozostaje otwarte

- to nie jest pelne emergentne wyprowadzenie spinora z jednego nadsolitonu,
- gamma matrices i spin bundle sa tutaj dolaczane jako kontrolowane rozszerzenie,
- nadal brak theorem-level closure dla kernel-origin fermionow.

### Ocena

- zgodnosc z ontologia: `umiarkowana`
- zgodnosc z admissible strict references: `wysoka`
- gotowosc wykonawcza: `wysoka`
- rygor recenzencki dzis: `partial but operationally useful`

## Rozstrzygniecie A5

Na obecnym stanie repo najlepszy podzial jest taki:

1. `Route primary`
   - `spinor-emergent from 3D topology`
   - bo to jest jedyna droga zgodna z ontologia jednego nadsolitonu,
   - ale pozostaje hipoteza konstrukcyjna, nie lokalnie dowiedziony pakiet.
2. `Route control / fallback`
   - `minimal spin-bundle extension`
   - bo daje szybszy i czystszy pomost do dalszych warstw `A6-A8`,
   - i opiera sie bardziej na admissible strict references niz na legacy heuristics.

To znaczy:
- glowny program konstrukcyjny idzie emergentnie,
- ale warstwa kontrolna pozostaje gotowa jako conservative bridge.

## Macierz statusu po A5

| Obiekt | Status | Uwagi |
|---|---|---|
| 3D topologiczny origin spinu `1/2` | `hypothesis branch with partial support` | legacy support istnieje, ale nie jest traktowany jako dowod |
| odrzucenie 2D vortex jako drogi do fermionow | `negative-control signal` | sygnal legacy, nie finalny proof artifact |
| Dirac-like action block | `partially derived / extension-backed` | `QW-2121/2127` |
| gamma matrices z nadsolitonowej liniaryzacji | `unresolved` | brak jawnego derivation package |
| anomaly/charge closure | `partially derived` | `QW-2189` |
| kernel-mode gauge scaffold | `partially derived` | `QW-2190` |
| pelna unikalnosc representation map | `blocked` | `QW-2191` |
| minimal spin-bundle extension | `available as control route` | nie emergentne ontologicznie |

## Co A5 rzeczywiscie zamyka

`A5` zamyka:
- jawny wybor dwoch dopuszczalnych drog spinorowych,
- degradacje badan legacy do roli heurystyki / negatywnej kontroli,
- odrzucenie 2D vortex route jako niewystarczajacej co najwyzej na poziomie sygnalu metodologicznego,
- uznanie 3D topological route za glowna droge zgodna z teoria,
- utrzymanie minimal spin-bundle extension jako kontrolnego mostu.

## Co pozostaje otwarte po A5

- jawny derivation package dla gamma matrices,
- jawny Dirac operator wyprowadzony z liniaryzacji `A1-A4`,
- pelna integracja spinor route z gauge reconstruction,
- pelna integracja spinor route z gravity bridge,
- theorem-level closure `L3/L18/L19`.

## Anti-overclaim

`A5` nie twierdzi, ze:
- spinory Diraca zostaly juz wyprowadzone,
- macierze gamma zostaly juz wyprowadzone,
- `SU(3)xSU(2)xU(1)` zostalo juz unikalnie zrekonstruowane z jednego nadsolitonu,
- anomaly package jest juz w pelni kernel-origin i axiom-free,
- most do GR jest juz domkniety.

## Produkt etapu

- jawny split `emergent route / control extension route`,
- matrix `derived / partially derived / blocked / control-only`,
- lista poprzednich badan w repo, ktore rzeczywiscie zasilaja dalszy tor.

## Nastepny krok

Naturalnym kolejnym ruchem jest `A6`:
- wziasc `Route primary` z `A5`,
- zestawic go z kernel-mode gauge scaffold `QW-2190`,
- i zbudowac jawna warstwe `gauge reconstruction` bez falszywego claimu unikalnosci.
