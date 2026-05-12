# P1315 — O3.2 R_common normalization and translation maps (strict-only, PL)

Status: `O3_2_EXECUTED_NORMALIZATION_SPEC_PARTIAL`
As of: `2026-05-12`
Depends on: `P1313`, `P1314`

## Cel
Wykonać O3.2 z `P1313` w sposób formalny i audytowalny:
1. zdefiniować wspólną reprezentację `R_common_v1`,
2. zdefiniować translacje `T_1..T_4` dla klas C1–C4,
3. zdefiniować test jednoznaczności translacji i status PASS/FAIL.

## Założenie rygoru
- Brak claimu closure `QW-2191`.
- Brak transferu legacy->strict.
- Brak mieszania strict-only z NB-only.

## Definicja reprezentacji `R_common_v1`
Każdy kandydat po translacji jest reprezentowany jako rekord:

```text
R_common_v1(c) := (
  branch_sign,            # {-1,+1}
  phase_class,            # [0,pi) modulo orientacyjne symetrie strict
  amplitude_class,        # klasa amplitudy znormalizowana do przedziału [0,1]
  selector_slot_residual, # {closed,open} + jawny tag źródła residualu
  admissibility_tag       # {admissible,rejected}
)
```

### Reguły kanoniczne
1. `branch_sign` ustalamy po orientacji strict i minimalnej regule tie-break:
   jeśli dwie orientacje są równoważne gauge-like, wybieramy reprezentanta o
   najmniejszej normie perturbacyjnej.
2. `phase_class` redukujemy modulo dozwolone symetrie strict (nie przez
   arbitralny gauge fix spoza klasy).
3. `amplitude_class` jest obiektem porównawczym (nie absolutnym źródłem
   ontologicznym).
4. `selector_slot_residual` musi jawnie wskazać, czy residual slot został
   zamknięty wewnętrznie strict, czy pozostaje otwarty.

## Mapy translacji `T_i: C_i -> R_common_v1`

### T1 dla C1 (projective selector-state)
- wejście: reprezentant stanu projektowego,
- wyjście:
  - `branch_sign`: znak orientacji po redukcji projektowej,
  - `phase_class`: klasa fazy projektowej,
  - `amplitude_class`: znormalizowana skala projektowa,
  - `selector_slot_residual`: zwykle `open` (chyba że jawnie domknięte strict),
  - `admissibility_tag`: `admissible`.

### T2 dla C2 (directed premise-based)
- wejście: directed selector-state,
- wyjście:
  - `branch_sign`: znak z directed branch signal,
  - `phase_class`: klasa fazy po normalizacji kierunkowej,
  - `amplitude_class`: skala directed po normalizacji,
  - `selector_slot_residual`: `open` lub `closed` zależnie od jawnego dowodu,
  - `admissibility_tag`: `admissible`.

### T3 dla C3 (sigma-int theta-pair ingredient)
- wejście: theta-pair ingredient + mode assignment,
- wyjście:
  - `branch_sign`: znak osi po translacji theta-pair,
  - `phase_class`: klasa fazowa ingredientu,
  - `amplitude_class`: skala ingredientowa,
  - `selector_slot_residual`: jawny tag `open(Z2/eps)` jeśli slot nieusunięty,
  - `admissibility_tag`: `admissible` (z residual flag).

### T4 dla C4 (observer-limit induced)
- wejście: observer-limit readout,
- wyjście:
  - `branch_sign`: znak kierunku readout,
  - `phase_class`: klasa fazy odczytu limitowego,
  - `amplitude_class`: znormalizowana odpowiedź readout,
  - `selector_slot_residual`: domyślnie `open`, chyba że strict proof mówi inaczej,
  - `admissibility_tag`: `admissible`.

## Dozwolone transformacje gauge-like (lista zamknięta v1)
`G_allowed_v1`:
1. redukcja fazy modulo symetrie strict,
2. równoważność skali amplitudowej przez dodatnie przeskalowanie porównawcze,
3. orientacyjne odwrócenie z kompensacją fazy, jeśli zachowuje ten sam stan
   fizyczny w sensie strict.

Transformacje spoza `G_allowed_v1` są niedozwolone w O3.2.

## Test jednoznaczności translacji (O3.2-UNIQUENESS)
Dla każdego `c in C_i` i dla dowolnych dwóch reprezentantów `r, r'` tej samej
klasy równoważności strict:

```text
normalize(T_i(r)) == normalize(T_i(r'))
```

w granicach tolerancji `tau_strict_v1`.

### Status logiczny
- `O3.2 = PASS`, gdy test jednoznaczności przechodzi dla C1–C4 i nie ma
  niejawnych transformacji spoza `G_allowed_v1`.
- `O3.2 = FAIL`, gdy istnieje choć jeden kandydat z niejednoznaczną translacją.

## Wynik na dziś
- `R_common_v1`: `DEFINED`
- `T_1..T_4`: `DEFINED`
- `G_allowed_v1`: `DEFINED_CLOSED_LIST`
- `O3.2-UNIQUENESS formal statement`: `DEFINED`
- `O3.2 execution replay`: `NOT_YET_RUN`

**Aktualny werdykt:** `O3.2_SPEC_READY_BUT_NUMERICAL_REPLAY_PENDING`.

## Konsekwencja dla QW-2191
Ten pakiet nadal **nie** domyka `QW-2191`. On tylko usuwa niejawność mapowania
przed testami parowymi O3.3 i counterexample sweep O3.4.

## Rekomendowany następny uczciwy krok
Wykonać **O3.3**: pairwise direction-separation na C1–C4 po `R_common_v1`, tj.
macierz porównań `(Ci,Cj)` z decyzją: `equivalent` albo `inadmissible branch`.

## Dla laika
Ustaliliśmy jedną wspólną „tabelkę wyników” dla wszystkich legalnych kompasów i
reguły tłumaczenia każdego kompasu do tej tabelki. Teraz można uczciwie porównać
kompasy para po parze i sprawdzić, czy naprawdę pokazują to samo.
