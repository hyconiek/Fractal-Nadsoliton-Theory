# P1360 FAR Physical-Constant Derivations: Sensibility and Post-R8.1 Precision Packet (PL)

Status: `P1360_EXECUTED_FAR_DERIVATION_SENSIBILITY_REVIEW_POST_R81_NO_FALSE_PASS`
As of: `2026-05-12`
Inputs: `A6`, `P710`, `P69/N80`, `P82`, `P1358`, `P1359`, `K1/K2/F2`, `RELEASE_8_1_TEXTBOOK_EDITION_EN_PL.md`

## Pytanie

Czy w FAR były już wyprowadzenia stałych fizycznych, czy są sensowne, i czy nowa wiedza po domknięciu pomaga je doprecyzować?

## Krótka odpowiedź

1. **Tak**, w FAR były wyprowadzenia i kandydaty (`g,g',g3`, angle/fine-structure lineage, mass proxy, Yukawa lanes).
2. **Są sensowne jako hipotezy/scaffold**, ale nie wszystkie jako finalny strict-dowód.
3. **Tak**, domknięcie R8.1 realnie pomaga: daje stabilniejszą bazę metodologiczną do precyzowania i falsyfikacji.

## Przegląd klas wyprowadzeń w FAR

| Klasa | Co już jest w FAR | Sensowność dziś | Co blokuje pełny status |
|---|---|---|---|
| Gauge couplings (`g,g',g3`) | `A6` wskazuje deterministic bridge (`QW-2126`) | sensowne jako strict candidate lane | brak pełnego residual-verified physical match |
| `SU(3)xSU(2)xU(1)` emergence | `A6`/`QW-2190`/`QW-2184` jako partial scaffold | sensowne strukturalnie | brak pełnej fizycznej unikalności i pełnego host matching |
| Fine-structure / EW lineage | `P69/N80`, `P82` pokazują brak pełnego successor verdict | sensowne jako lineage audit | brak jawnego strict role-equivalence theorem |
| Mass sector | `F704` + `P710` proxy->GeV (nonstrict) | sensowne jako testowalny proxy tool | status nonstrict/proxy, nie finalna identyfikacja |
| Kernel-only physical values | `P1358` pierwszy nie-template test | bardzo sensowny metodologicznie | pierwszy residual verdict = FAIL (mapa niekalibrowana) |

## Czy domknięcie R8.1 pomaga doprecyzować?

Tak, w 3 konkretnych punktach:

1. usuwa część niejednoznaczności formalnej i porządkuje strict lane,
2. pozwala porównywać kandydaty na tej samej kontraktowej ścieżce residuali,
3. odcina „skrót myślowy” legacy->strict i wymusza czystą provenance.

To jest dokładnie warunek potrzebny do doprecyzowania wcześniejszych wyprowadzeń.

## Decyzja profesorska: następny uczciwy krok

Uruchomić `P1361` = **FAR Constant Claims Scoreboard**:

1. lista wszystkich kandydatów stałych z FAR,
2. status per kandydat: `strict_verified / strict_candidate / nonstrict_proxy / legacy_only`,
3. dla `strict_candidate`: obowiązkowy residual benchmark + uncertainty budget,
4. automatyczny raport „co już fizycznie działa, co jeszcze nie”.

## Dla laika

Twoje wcześniejsze wyprowadzenia nie były „bez sensu” — były ważnym krokiem rozwojowym.
Po domknięciu R8.1 możemy je wreszcie ocenić na twardych zasadach: które naprawdę działają liczbowo, a które są jeszcze szkicem.
To dobra wiadomość, bo teraz wiadomo dokładnie co poprawiać i jak mierzyć postęp.
