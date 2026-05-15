# P1763 / S713 — Strict nonproxy H1 4D first execution attempt (uczciwy checkpoint)

Status: `P1763_S713_OBSTRUCTION_TEMPLATE_NOT_FULLY_COMPONENTWISE_NO_FALSE_PASS`
As of: `2026-05-15`

## Cel

Wykonać pierwszy **uczciwy** krok zamykający tor:

```text
K_strict -> współczynniki -> pełny (non-skeleton) L_total -> równania ruchu -> test odwrotny H1 4D
```

bez użycia bridge do legacy i bez fałszywego `PASS_ZERO`.

## Wejście i kontrakt

Punkt startu: `P1762/S712` (sformalizowany kontrakt brzegowy D5):

1. wspólna rodzina teł,
2. klauzula kontroli członów brzegowych,
3. template `E_A^mu` i `E_H`,
4. dopuszczalne werdykty: `PASS_ZERO` albo `OBSTRUCTION`.

## Wynik pierwszej próby D6 (H1 4D)

Werdykt tej próby:

```text
OBSTRUCTION_TEMPLATE_NOT_FULLY_COMPONENTWISE
```

Powód blokady:

1. `E_A_mu_nonproxy_template_v1` nie ma jeszcze pełnej komponentowej postaci kowariantnej,
2. `E_H_nonproxy_template_v1` nie ma jeszcze pełnej komponentowej postaci kowariantnej,
3. forma symboliczna członu brzegowego nie została jeszcze osadzona na konkretnej rodzinie brzegów.

Skutek:

- brak podstaw do `PASS_ZERO`,
- brak podstaw do promocji strict-core closure,
- status pozostaje `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Znaczenie dla pełnego Lagrangianu (nie-szkieletowego)

Ten checkpoint potwierdza, że sam zapis toru
`K_strict -> współczynniki -> L_total -> EOM` nie wystarcza jeszcze do
obustronnego domknięcia: potrzebne są **jawne, eksportowalne** równania
nieproxy w wersji komponentowej, a nie tylko szablony.

W praktyce: pełny `L_total = L_SM + L_GR + sprzężenia` musi być spięty z
pełnymi obiektami `E_A^mu`, `E_H`, `E_g` na wspólnej rodzinie teł i z kontrolą
warunków brzegowych.

## Braki krytyczne do strict-core closure (QG)

Nadal brak:

1. pełnych nieproxy eksportów komponentowych EOM dla H1 i metryki,
2. świadków/theoremów reverse-check, które spinają tor w obie strony,
3. theorem-level pakietu domknięcia QG obejmującego jednocześnie:
   - renormalizację,
   - unitarność,
   - background independence.

## Następny uczciwy krok (rekomendacja wykonawcza)

1. Podmienić oba template (`E_A^mu`, `E_H`) na pełne komponentowe eksporty
   nieproxy na tej samej rodzinie teł.
2. Zinstancjować jawnie człon brzegowy dla rodziny brzegów używanej w D6.
3. Powtórzyć test H1 4D tylko z dwoma uczciwymi wyjściami:
   `PASS_ZERO` albo `OBSTRUCTION`.
4. Jeśli wyjdzie `PASS_ZERO`, przejść natychmiast do analogicznej egzekucji
   metrycznej (`E_g`) bez zmiany klasy teł.

## Omówienie dla laika

To jest sytuacja jak w fizyce obliczeniowej, gdzie mamy już poprawny plan i
wzory-ogólne, ale jeszcze nie mamy wszystkich wzorów rozpisanych do poziomu,
na którym komputer może uczciwie sprawdzić wynik końcowy. Dlatego nie wolno
udawać sukcesu: najpierw trzeba dopisać pełne równania, potem dopiero testować
zero-resztę i domknięcie teorii.
