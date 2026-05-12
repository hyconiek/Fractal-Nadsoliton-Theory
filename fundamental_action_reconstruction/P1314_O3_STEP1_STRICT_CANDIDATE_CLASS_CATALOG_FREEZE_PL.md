# P1314 — O3.1 Strict candidate-class catalog freeze for QW-2191 (PL)

Status: `O3_1_EXECUTED_CATALOG_FROZEN_PARTIAL`
As of: `2026-05-12`
Depends on: `P1312`, `P1313`

## Cel
Wykonać pierwszy realny krok z `P1313` (`O3.1`) w sposób audytowalny:
- zamrozić skończony katalog klas kandydatów strict-only,
- jawnie oznaczyć interfejs, constraints i status dopuszczalności,
- nie ogłaszać jeszcze `QW-2191 = CLOSED`.

## Zakres i rygor
Ten pakiet klasyfikuje klasy kandydatów do testu unikalności kierunku w strict lane.
Nie wykonuje jeszcze O3.2–O3.5.

## Katalog klas `C_strict` (zamrożony na dziś)

### C1 — Projective selector-state class
**Opis:** klasa projektowa (projective) oparta o globalne stany selectorowe bez pełnego upgrade do strict-core closure.

**Interfejs roboczy:**
- wejście: lokalne/kanoniczne obiekty orientacji i mapowania projektowego,
- wyjście: kierunek projektowy modulo symetrie projektowe.

**Constraints:**
- brak claimu strict-core closure,
- brak claimu globalnego discharge `QW-2191`.

**Status admissibility:** `ADMISSIBLE_AS_CANDIDATE_CLASS`.

### C2 — Directed premise-based selector-state class
**Opis:** klasa kierunkowa (directed) z jawną premią symetryzacji łamanej przez premise-based directed state.

**Interfejs roboczy:**
- wejście: premise-based directed selector state,
- wyjście: kierunek zorientowany (directed branch signal).

**Constraints:**
- nie jest automatycznie strict-core closure,
- wymaga jawnego oddzielenia od NB i od rozszerzeń aksjomatycznych.

**Status admissibility:** `ADMISSIBLE_AS_CANDIDATE_CLASS`.

### C3 — Sigma-int theta-pair ingredient class
**Opis:** klasa ingredientowa `sigma-int -> theta-pair` z residualnym slotem `Z2/eps` w granicach strict.

**Interfejs roboczy:**
- wejście: strict ingredient (`theta pair` / mode assignment),
- wyjście: kandydat kierunku lokalnego do transportu na wspólną reprezentację.

**Constraints:**
- brak automatycznego zamknięcia slotu selektora,
- brak claimu globalnej unikalności.

**Status admissibility:** `ADMISSIBLE_WITH_RESIDUAL_SLOT_FLAG`.

### C4 — Observer-limit readout induced class
**Opis:** klasa kierunków indukowanych przez strict observer-limit readout operator.

**Interfejs roboczy:**
- wejście: strict-core readout operator,
- wyjście: kierunek odczytu na warstwie limitowej.

**Constraints:**
- observer-limit operator nie implikuje samodzielnie selector closure,
- brak claimu discharge `QW-2191`.

**Status admissibility:** `ADMISSIBLE_AS_CANDIDATE_CLASS`.

### C5 — Non-admissible closure-claim class (wykluczona)
**Opis:** wszystkie klasy, które deklarują `QW-2191` closure przez samo mapowanie formalne bez rozładowania przeszkody unikalności.

**Status admissibility:** `REJECTED_NONADMISSIBLE`.

## Zamrożona tabela O3.1
| Klasa | Kod | Interfejs jawny | Residual slot | Dopuszczalność |
|---|---|---|---|---|
| Projective selector-state | C1 | tak | możliwy | admissible |
| Directed premise-based | C2 | tak | możliwy | admissible |
| Sigma-int theta-pair ingredient | C3 | tak | jawny `Z2/eps` | admissible with flag |
| Observer-limit induced | C4 | tak | możliwy | admissible |
| Closure-claim shortcut | C5 | n/a | n/a | rejected |

## Decyzja O3.1
- `O3.1 = PASS (catalog freeze)` dla zakresu klasyfikacyjnego.
- `O3 global != PASS` (O3.2–O3.5 pozostają niewykonane).
- `nonuniqueness_residual` pozostaje `OPEN`.

## Warunki przejścia do O3.2
Aby uruchomić O3.2 (normalizacja reprezentacji), trzeba dla C1–C4 dopisać:
1. jednolity format reprezentacji `R_common_v1`,
2. funkcję translacji `T_i: C_i -> R_common_v1`,
3. listę dozwolonych transformacji gauge-like do zmodowania.

## Czego ten dokument nie twierdzi
- Nie twierdzi domknięcia `QW-2191`.
- Nie twierdzi strict-core selector closure.
- Nie twierdzi ToE closure.

## Rekomendowany następny uczciwy krok
Wykonać **O3.2**: sformalizować `R_common_v1` i jawne translacje `T_1..T_4`, a następnie udowodnić ich jednoznaczność w tolerancji strict.

## Dla laika
Zrobiliśmy listę legalnych typów „kompasów” i odrzuciliśmy fałszywy skrót "ogłaszam domknięcie bez dowodu". Następny etap to przetłumaczyć wszystkie legalne kompasy na jedną wspólną „skalę”, żeby dało się uczciwie porównać wskazania.
