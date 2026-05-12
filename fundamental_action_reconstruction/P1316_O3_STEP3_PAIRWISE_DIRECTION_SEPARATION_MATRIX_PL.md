# P1316 — O3.3 Pairwise direction-separation matrix (strict-only, PL)

Status: `O3_3_EXECUTED_PAIRWISE_MATRIX_PARTIAL`
As of: `2026-05-12`
Depends on: `P1314`, `P1315`

## Cel
Wykonać O3.3: porównać parowo klasy C1–C4 po translacji do `R_common_v1` i
ustalić, czy:
1. są operacyjnie równoważne kierunkowo (`equivalent`), albo
2. jedna gałąź jest niedopuszczalna (`inadmissible branch`), albo
3. pozostaje nierozstrzygnięta rozbieżność wymagająca O3.4.

## Wejście formalne
- Katalog klas C1–C4 z `P1314`.
- Reprezentacja `R_common_v1` i mapy `T_1..T_4` z `P1315`.

## Reguła porównania pary `(Ci, Cj)`
Dla każdej pary liczymy relację:

```text
Delta_ij := compare(normalize(T_i(Ci)), normalize(T_j(Cj)))
```

z wynikiem w zbiorze:
- `EQUIVALENT_DIRECTION`,
- `INADMISSIBLE_BRANCH`,
- `RESIDUAL_AMBIGUITY`.

## Macierz O3.3 (stan bieżący)

| Para | Werdykt | Uzasadnienie strict-only |
|---|---|---|
| (C1, C2) | `RESIDUAL_AMBIGUITY` | projective vs directed daje zgodny sygnał orientacyjny, ale bez pełnego dowodu eliminacji slotu residualnego |
| (C1, C3) | `RESIDUAL_AMBIGUITY` | C3 posiada jawny residual `Z2/eps`, którego nie można jeszcze zredukować do jednoznacznego kierunku globalnego |
| (C1, C4) | `RESIDUAL_AMBIGUITY` | observer-limit wspiera kierunek, ale nie domyka selector-source uniqueness |
| (C2, C3) | `RESIDUAL_AMBIGUITY` | directed branch i theta-ingredient są zgodne lokalnie, brak globalnej eliminacji kontr-gałęzi |
| (C2, C4) | `RESIDUAL_AMBIGUITY` | zgodność sygnału kierunku nie jest równoznaczna z domknięciem przeszkody QW-2191 |
| (C3, C4) | `RESIDUAL_AMBIGUITY` | residual slot w C3 pozostawia alternatywną gałąź dopuszczalną do czasu O3.4 |

## Wniosek O3.3
- Brak pary z werdyktem `INADMISSIBLE_BRANCH` na podstawie obecnych danych.
- Brak pełnej sieci `EQUIVALENT_DIRECTION` wymaganej do domknięcia O3.
- Wszystkie pary utrzymują `RESIDUAL_AMBIGUITY`.

**Status:** `O3.3 = PARTIAL_FAIL_FOR_CLOSURE`.

## Konsekwencja logiczna
- `nonuniqueness_residual` pozostaje `OPEN`.
- `QW-2191` pozostaje `NOT_CLOSED` w strict-core.
- Nie wolno ogłaszać strict closure na podstawie obecnej macierzy.

## Minimalny pakiet do O3.4 (counterexample sweep)
Aby przejść dalej uczciwie, trzeba uruchomić kontrprzykłady co najmniej w osiach:
1. perturbacje fazowe w granicach `tau_strict_v1`,
2. perturbacje amplitudowo-skalujące w `G_allowed_v1`,
3. perturbacje slotu residualnego `Z2/eps` (C3),
4. próby wymuszenia alternatywnego `branch_sign` bez naruszenia admissibility.

## Czego dokument nie twierdzi
- Nie twierdzi domknięcia O3.
- Nie twierdzi domknięcia `QW-2191`.
- Nie twierdzi globalnego ToE closure.

## Rekomendowany następny uczciwy krok
Wykonać **O3.4**: formalny counterexample sweep na wszystkich parach z
`RESIDUAL_AMBIGUITY`; celem jest wykazać, że alternatywne gałęzie są albo
niedopuszczalne, albo redukowalne do tej samej klasy kierunku.

## Dla laika
Porównaliśmy wszystkie legalne kompasy parami. Na razie żaden nie został
obalony, ale też nie mamy dowodu, że wszystkie są w 100% tym samym kompasem.
Dlatego trzeba teraz aktywnie szukać kontrprzykładów i sprawdzić, czy któryś
kompas naprawdę może wskazać inny legalny kierunek.
