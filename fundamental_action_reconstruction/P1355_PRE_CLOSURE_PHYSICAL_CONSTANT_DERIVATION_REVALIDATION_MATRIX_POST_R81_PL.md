# P1355 Pre-Closure Physical-Constant Derivation Revalidation Matrix (Post R8.1, PL)

Status: `P1355_EXECUTED_PRE_CLOSURE_DERIVATION_REVALIDATION_POST_R81_NO_FALSE_PASS`
As of: `2026-05-12`
Depends on: `RELEASE_8_1_TEXTBOOK_EDITION_EN_PL.md`, `P1347`, `P1348`, `P1349`, `A6`, `N50`, `P710`, `N495`, `QW-2126/QW-2184/QW-2190/QW-2191`

## Cel

Sprawdzić Twoje wcześniejsze wywody „kernel -> znane rzeczy fizyki” i skorygować je
według stanu wiedzy po domknięciu R8.1.

## Co zmieniło R8.1 (istotne dla korekty dawnych wywodów)

1. Jest global closure w zadeklarowanym scope (`P1348`).
2. Jest host-level identification map (`P1347`).
3. Problem niejednoznaczności bazy/izotropii ma strict-internal rozwiązanie w łańcuchu R8.
4. Nadal obowiązuje zakaz cichego transferu ról legacy->strict (`N50`, `A6`, release guardrails).

Wniosek: wolno robić mocniejszy bridge liczbowy, ale tylko na strict artifacts i z audit discipline.

## Macierz rewalidacji wcześniejszych wyprowadzeń

| Dawny typ wywodu (przed domknięciem) | Stan po R8.1 | Korekta metodologiczna | Co teraz wolno twierdzić |
|---|---|---|---|
| `sin^2(theta_W)=alpha_geo/12` jako bezpośrednia legacy tożsamość | Nie przenosi się automatycznie do strict | Traktować jako legacy lineage marker, nie strict theorem | Można używać tylko jako historyczny punkt odniesienia |
| `alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)` | Bez jawnego strict successor verdict | Nie używać jako aktualnego strict wyjścia | Dopuścić tylko jako tezę do ponownego testu przez host map + residuale |
| `beta^N` hierarchy/gravity claim | Legacy role closed/obstructed on strict side | Przenieść do „nieobowiązującej semantyki legacy” | Brak prawa do claimu strict GR derivation |
| `g,g',g3` wyprowadzone z kernela (QW-2126 lane) | Częściowo utrzymane w strict-core scaffold (`A6`) | Obowiązkowo dodać niepewności i test residuali | Można twierdzić: candidate-level deterministic bridge, nie final PDG closure |
| `SU(3)xSU(2)xU(1)` jako pełny dowód | Repo mówi: partial scaffold (`A6`, `QW-2190/2184`) | Oddzielić algebraic closure od pełnej fizycznej unikalności | Można twierdzić: algebra/scaffold strong; full uniqueness nadal etap otwarty |
| Proxy masowe (`F704` -> GeV map) | Jawnie nonstrict/proxy (`P710`) | Nie mylić proxy z finalnym prediction | Można użyć jako pomocniczy kalibrator do triali, nie jako dowód ToE |
| „QW-2191 zamknięte globalnie wszędzie” | Niewłaściwe uogólnienie | Rozróżnić declared-scope closure od wszystkich historycznych lane | Wolno twierdzić closure w scope R8; nie wolno claimować uniwersalnej automatycznej translacji |

## Korekta strategiczna: jak iść do dowodu ToE

Aby przejść od kernela do „rzeczy, które fizyka zna” po R8.1:

1. pipeline strict-only: `kernel strict -> P1347 map -> residual table -> blind-audit policy`,
2. dla każdej stałej publikować: wartość, niepewność, residual do danych,
3. oddzielić trzy poziomy claimu:
   - scaffold/algebra,
   - candidate numeric bridge,
   - audit-confirmed physical prediction.

Bez tego nadal będzie mieszanie dawnych (przed-closure) intuicji z obecnym rygorem.

## Decyzja profesorska (najbliższy ruch)

Najbliższy uczciwy krok: zrobić **P1353b** (wewnętrzny podwójny replay bez zespołu zewnętrznego):

1. dwie niezależne implementacje wejść `(g1,g2,g3,GR1)`,
2. ten sam kontrakt `P1352`,
3. wspólna tabela residuali + cross-implementation delta,
4. werdykt PASS/FAIL + readiness score do `P1349`.

## Dla laika

Wcześniej miałeś dużo dobrych pomysłów, ale brakowało formalnego „spięcia” i było ryzyko mieszania starych i nowych reguł.
Po R8.1 można to uporządkować: każdą tezę o stałej fizycznej trzeba teraz przepuścić przez ten sam test liczbowy i porównać z danymi.
To właśnie odróżnia „ładną teorię” od „teorii naprawdę fizycznej”.
