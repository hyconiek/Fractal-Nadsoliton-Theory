# A3 Kernel Analysis

Status: `A3_EXECUTED_MINIMAL_BRANCH_OPERATOR_SPLIT`
As of: `2026-03-06`

## Cel

Wykonac pierwszy uczciwy etap analizy kernela na dokladnie tej samej galezi, ktora zostala wykonana w `A2`:
- `single-foundation`,
- `gauge-off`,
- `metric-spectator`.

To nie jest jeszcze pelna analiza calej teorii. To jest analiza operatora drugiej wariacji na minimalnej galezi radialnej, ktora ma stanowic poprawne wejscie do `A4`.

## Wejscie z A2

`A2` dostarczylo:
- tlo `Psi^A(x) = rho(r) n^A(Omega)`,
- `Phi(x) = phi(r)`,
- zredukowana akcje radialna,
- matchingowe rownania E-L,
- warunki brzegowe i klasyfikacje `forced / optional / gauge-choice-dependent`.

`A3` bierze dokladnie te dane i nie rozszerza jeszcze modelu o:
- aktywny sektor gauge,
- backreaction grawitacyjny,
- fermiony.

## Fluktuacje wokol tla

Rozwijamy tlo jako:

```text
rho(r) -> rho_bg(r) + delta rho(r)
phi(r) -> phi_bg(r) + delta phi(r)
n^A(Omega) -> n_bg^A(Omega) + delta n_perp^A(Omega)
```

Na wykonanej galezi:
- `A_mu^I = 0`,
- `g_{mu nu} = eta_{mu nu}`,
- wiec podstawowy blok fizyczny zaczyna sie od dubletu:
  - `xi = (delta rho, delta phi)^T`.
  - `delta rho` jest bezposrednia fluktuacja pola fundamentalnego,
  - `delta phi` jest efektywna fluktuacja warstwy porzadku zwiazanej z ta sama ontologia jednego nadsolitonu.

## Jawny operator drugiej wariacji

W minimalnej radialnej redukcji operator fizyczny ma postac macierzowego operatora Sturm-Liouville:

```text
O_phys = - d/dr [ K_2(r) d/dr ] + M_2(r)
```

z:

```text
K_2(r) = r^(d-1) [[ K(rho_bg,phi_bg), 0 ],
                  [ 0,                 1 ]]
```

oraz:

```text
M_2(r) =
  r^(d-1) Hess_(rho,phi) V_eff |_(bg)
  - r^(d-3) Hess_(rho,phi) T_top |_(bg)
  + connection/background-gradient corrections from K(rho,phi).
```

Wniosek operacyjny:
- fizyczny kernel nie jest skalarna liczba,
- jest macierzowym operatorem z mieszaniem kanalow `delta rho` i `delta phi`,
- dlatego jakikolwiek prosty claim typu "kernel jest dodatnio okreslony" bez projekcji i rozdzialu modow bylby zbyt gruby.

## Obowiazkowy rozdzial modow

### 1. Zero modes

Na tej galezi oczekiwane sa:
- translacyjne zero modes, jesli rzeczywiscie istnieje zlokalizowane rozwiazanie,
- orientacyjne zero modes, jesli `n^A` ma ciagly manifold modulow,
- scale mode tylko wtedy, gdy tla sa na galezi krytycznej/BPS-like; to nie jest tutaj zakladane.

### 2. Gauge modes

Na wykonanej galezi `gauge-off`:
- w aktualnym zredukowanym zbiorze zmiennych nie ma aktywnego fizycznego sektora gauge,
- wiec nie ma jeszcze pelnego pakietu gauge modes do odjecia.

Ale rygor wymaga dopisania granicy:
- po aktywacji `A_mu^I` pojawia sie czesc pure-gauge,
- i dopiero po jej projekcji wolno mowic o coercivity calego bosonicznego kernela.

### 3. Physical modes

Na minimalnej galezi fizyczne mody to:
- radialny kanal amplitudowy `delta rho`,
- kanal order-parameter `delta phi`,
- ich mieszanie przez `M_2(r)`,
- topologiczno-ksztaltowe fluktuacje ortogonalne do zero modes.

## Stabilnosc: poprawne kryteria sektorowe

### Bosonic sector (Euclidean)

Na wykonanej galezi kryterium ma postac:
- najpierw odjac zero modes,
- nastepnie odjac gauge modes, jesli sektor gauge zostanie aktywowany,
- dopiero potem sprawdzac coercivity / positivity operatora `O_phys`,
- oraz brak tachyonic instability w zadeklarowanym scope.

### Fermionic sector

`A3` nie zamyka sektora fermionowego. Pozostaje poprawna granica:
- nie wolno zadac "dodatniosci operatora Diraca",
- trzeba docelowo sprawdzac self-adjointness / ellipticity / positivity of `D^dagger D` / reflection positivity.

### Lorentzian sector

`A3` nie przechodzi jeszcze do pelnej sygnatury lorentzowskiej. Poprawna granica pozostaje:
- well-posedness,
- hyperbolicity,
- bounded-below Hamiltonian,
- brak ghostow.

## Co A3 rzeczywiscie ustala

`A3` ustala:
- jawna forme fizycznego operatora drugiej wariacji na minimalnej galezi,
- to, ze operator jest blokiem macierzowym, a nie skalarna liczba,
- oczekiwany split `zero / gauge / physical modes`,
- wymaganie, ze kazdy przyszly claim o stabilnosci musi byc liczony po projekcjach.

## Co pozostaje otwarte po A3

- czy `O_phys` jest coercive dla konkretnego tla i konkretnego wyboru `K`, `V_eff`, `T_top`,
- czy po wlaczeniu sektora gauge pojawiaja sie dodatkowe niestabilnosci,
- czy orientacyjne mody `delta n_perp^A` sa modami zerowymi, czy niosa dodatnia mase efektywna,
- jak wyglada kernel po wlaczeniu fermionow i grawitacji,
- czy analiza po coarse-graining prowadzi do rzeczywistego emergentnego RG.

## Anti-overclaim

`A3` nie twierdzi, ze:
- globalna stabilnosc jest juz udowodniona,
- fermionowy kernel jest juz zamkniety,
- sektor gauge jest juz rozladowany,
- GR-limit jest juz kontrolowany,
- jest juz theorem-level positivity/well-posedness package.

## Produkt etapu

- jawny operator drugiej wariacji `O_phys`,
- macierz `K_2(r)` i efektywny blok masowy `M_2(r)`,
- jawny split `zero / gauge / physical modes`,
- lista projection-before-claim constraints dla przyszlych etapow.

## Nastepny krok

Naturalnym kolejnym ruchem jest `A4`:
- zdefiniowac coarse-graining dla operatora `O_phys`,
- sprawdzic, ktore sprzezenia sa `emergent / inserted by hand / unresolved`,
- i dopiero wtedy wracac do twierdzen o RG emergence.
