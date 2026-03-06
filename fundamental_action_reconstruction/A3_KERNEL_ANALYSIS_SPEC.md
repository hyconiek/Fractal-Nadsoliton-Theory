# A3 Kernel Analysis Spec

Status: `READY_NOT_EXECUTED`

## Cel

Policzyc druga wariacje dzialania wokol tla supersolitonowego i zidentyfikowac operator fluktuacji.

## Obowiazkowy rozdzial modow

1. zero modes,
2. gauge modes,
3. physical modes.

Brak tego rozdzialu oznacza, ze dowolny claim o "dodatniosci kernela" jest fizycznie zbyt gruby.

## Poprawne kryteria sektorowe

### Bosonic sector (Euclidean)
- coercivity / positivity po gauge-fixingu,
- usuniecie zerowych modow,
- brak tachyonic instability w zadeklarowanym scope.

### Fermionic sector
- nie zadac "dodatniosci operatora Diraca",
- sprawdzac self-adjointness / ellipticity / positivity of `D^dagger D` / reflection positivity.

### Lorentzian sector
- well-posedness,
- hyperbolicity,
- bounded-below Hamiltonian,
- brak ghostow.

## Produkt etapu

- jawny operator drugiej wariacji,
- blokowa dekompozycja sektorow,
- tabela wymagan stabilnosci,
- lista modow wymagajacych gauge-fixingu lub constraints.
