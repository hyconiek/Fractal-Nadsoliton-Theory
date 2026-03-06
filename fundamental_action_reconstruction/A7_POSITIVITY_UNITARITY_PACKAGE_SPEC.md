# A7 Positivity / Unitarity Package

Status: `A7_EXECUTED_STRICT_SCOPE_PARTIAL_WITH_TERMINAL_POSITIVITY_AND_UNITARITY_OBLIGATIONS_OPEN`
As of: `2026-03-06`

## Cel

Po `A1-A6` zbudowac uczciwy pakiet `positivity / unitarity`, ale tylko na tym, co repo ma juz:
- jako strict-admissible internal evidence,
- jako branch-scope albo strict-scope closure,
- bez wprowadzania falszywego theorem-level PASS.

`A7` ma odpowiedziec nie na pytanie:
- "czy pelna globalna kwantyzacja, unitarity i reconstruction sa juz domkniete?"

tylko na pytanie:
- "jaki najmocniejszy pakiet positivity / causality / partial unitarity wolno juz utrzymac po `A1-A6`, z jawnym rozdzialem local-vs-global?"

## Strict-admissible warstwa wejscia

`A7` uzywa tylko:
1. `A3`
   - poprawne kryteria sektorowe dla kernela:
   - bosonic Euclidean: coercivity po projekcji `zero / gauge`,
   - fermionic: self-adjointness / ellipticity / positivity of `D^dagger D` / reflection positivity,
   - Lorentzian: well-posedness / hyperbolicity / bounded-below Hamiltonian / brak ghostow.
2. `QW-2186`
   - strict branch-scope spectral stability margin dla `A = K_total + m0^2 I`.
3. `QW-2133`
   - free-sector microcausality closure in scope.
4. `QW-2134`
   - perturbative interacting microcausality conditions.
5. `QW-2138`
   - proof-completion gate dla interacting microcausality z jawna granica all-orders machine-check.
6. `QW-2202`
   - strict local action + microcausality + renormalization schema zintegrowane; globalny pakiet QFT pozostaje otwarty.
7. `QW-2214`
   - terminal constructive positivity-to-reconstruction obligation `L5_O1a_O1`.
8. `QW-2216`
   - terminal unitary-scattering obligation `L5_O1b_O1`.

## Jawnie wykluczone z rdzenia A7

`A7` nie liczy jako zamkniety proof input:
- legacy / exploratory prior art,
- theorem-provider chain dla pelnego `L5` po warstwie globalnej de-axiomatization,
- `A8+` gravity / full SM+GR bridge,
- jakiegokolwiek domniemania, ze fermionic reflection-positivity i Lorentzian unitarity sa juz automatycznie rozladowane przez sam lokalny stack.

## Co A7 rzeczywiscie ustala

### 1. Bosonic Euclidean positivity

Z `A3` i `QW-2186` wolno utrzymac:
- poprawne kryterium dodatniosci nie jest skalarne i nie wolno go redukowac do sloganu "kernel jest dodatni",
- najpierw trzeba odjac `zero modes`,
- po aktywacji sektora gauge trzeba odjac `gauge modes`,
- dopiero potem wolno badac coercivity / positivity operatora bosonicznego,
- dla `A = K_total + m0^2 I` istnieje strict branch-scope certyfikat Weyla:
  - `lambda_min(A)=0.331677209978`,
  - perturbacje symetryczne z `||Delta||_2 < lambda_min(A)` zachowuja PSD.

To daje:
- realny branch-scope positivity margin,
- ale nie globalna coercivity dla nieograniczonych / nielokalnych / nieliniowych klas perturbacji.

### 2. Free i perturbative microcausality

Z `QW-2133` i `QW-2134` wolno utrzymac:
- free-sector microcausality przechodzi in scope,
- perturbative interacting microcausality conditions sa spelnione,
- brak pelnego all-orders constructive proof pozostaje jawny.

To daje:
- strict-scope local causality package,
- ale nie globalna all-orders closure.

### 3. All-orders microcausality scaffold

Z `QW-2138` wolno utrzymac:
- macierz `8/8` jawnych obligacji jest spelniona,
- jest jawna kontrola remainderu wysokich rzedow,
- nadal brak pelnego machine-checked distribution-level all-orders proof export.

To daje:
- silny proof-completion scaffold,
- ale nie theorem-level full constructive closure.

### 4. Strict local QFT stack

Z `QW-2202` wolno utrzymac:
- strict scope for local action + microcausality + renormalization schema is integrated and explicit,
- globalny pakiet QFT pozostaje jawnie otwarty w trzech komponentach:
  - nonperturbative existence/uniqueness,
  - global S-matrix unitarity,
  - reflection positivity / Wightman reconstruction.

To daje:
- uczciwy status `strict local stack integrated`,
- bez przejscia do global theorem-level closure.

### 5. Positivity-to-reconstruction boundary

Z `QW-2214` wolno utrzymac:
- `L5_O1a` jest zredukowane do jednego terminalnego twierdzenia:
  - `L5_O1a_O1`
- tresc obligacji:
  - constructive nonperturbative existence/uniqueness
  - z positivity-to-reconstruction theorem
  - dla complete FIN action.

To daje:
- jawny pojedynczy terminalny blocker dla warstwy `positivity -> reconstruction`,
- ale nie jego discharge.

### 6. Unitary scattering boundary

Z `QW-2216` wolno utrzymac:
- `L5_O1b` jest zredukowane do jednego terminalnego twierdzenia:
  - `L5_O1b_O1`
- tresc obligacji:
  - unitary S-matrix
  - asymptotic/scattering completeness
  - dla reconstructed nonperturbative QFT z complete FIN action.

To daje:
- jawny pojedynczy terminalny blocker dla warstwy unitarity/scattering,
- ale nie jego discharge.

## Macierz statusu A7

| Obiekt | Status | Zakres |
|---|---|---|
| bosonic Euclidean positivity margin | `partial_strict_branch_scope` | `QW-2186`, bounded symmetric perturbations |
| projection-before-claim kernel discipline | `established` | `A3` |
| free-sector microcausality | `partial_strict_scope` | `QW-2133` |
| perturbative interacting microcausality | `partial_conditional` | `QW-2134` |
| all-orders microcausality scaffold | `partial_proof_completion` | `QW-2138` |
| local action + microcausality + renorm stack | `integrated_strict_scope` | `QW-2202` |
| positivity-to-reconstruction theorem | `blocked_by_single_terminal_obligation` | `L5_O1a_O1`, `QW-2214` |
| unitary scattering completeness theorem | `blocked_by_single_terminal_obligation` | `L5_O1b_O1`, `QW-2216` |
| global reflection positivity / Wightman reconstruction | `open` | `QW-2202` |
| full Lorentzian ghost-free / hyperbolic closure | `unresolved_in_A7` | `A3` gives criteria only |

## Rozstrzygniecie A7

Najuczciwszy wynik `A7` brzmi:

- teoria idzie do przodu do poziomu **integrated strict-scope positivity/unitarity package**,
- istnieje realny branch-scope bosonic positivity margin,
- istnieje realny strict-scope causality stack,
- globalne theorem-level closure nadal nie istnieje,
- i pozostaje jawnie rozdzielone na co najmniej dwa terminalne blockery:
  - `L5_O1a_O1`
  - `L5_O1b_O1`

## Anti-overclaim

`A7` nie twierdzi, ze:
- globalna nonperturbative existence/uniqueness jest udowodniona,
- globalna reflection positivity / Wightman reconstruction jest udowodniona,
- globalny unitary S-matrix jest udowodniony,
- fermionowy sektor jest juz theorem-level domkniety,
- Lorentzowski pakiet hyperbolicity / ghost-freedom jest juz zamkniety,
- `L5` jest rozladowane.

## Produkt etapu

- jawna warstwa `positivity / unitarity package`,
- rozdzial:
  - `branch-scope certified`,
  - `strict-scope integrated`,
  - `terminal-theorem blockers`,
  - `still unresolved`.

## Nastepny krok

Naturalnym kolejnym ruchem jest `A8`:
- przejsc do `gravity bridge`,
- ale bez przenoszenia do tego etapu jakiegokolwiek falszywego claimu, ze `A7` juz rozladowalo globalny pakiet `L5`.
