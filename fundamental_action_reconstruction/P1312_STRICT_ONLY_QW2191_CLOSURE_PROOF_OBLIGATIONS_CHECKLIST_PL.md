# P1312 — Strict-only QW-2191 closure proof-obligations checklist (PL)

## Cel
Ten dokument definiuje **minimalny uczciwy pakiet dowodowy strict-only**, który musi zostać spełniony, aby:
1. `nonuniqueness_residual` mogło przejść na `CLOSED`,
2. `QW-2191` mogło przejść na `CLOSED` w lane strict-core,
3. bez mieszania tego z NB-only lub mostem legacy->strict.

## Założenia rygoru
- Brak cichego transferu ról z `K_legacy_ont` do `K_strict_gate`.
- Brak claimów globalnej closure poza strict-only lane.
- Każdy krok musi być reprodukowalny i audytowalny.

## Checklist proof-obligations (strict-only)

### O1. Selector-source existence (strict)
- [ ] Istnieje formalna teza o **istnieniu** wewnętrznego źródła selektora dla strict lane.
- [ ] Teza nie odwołuje się do axiom-tagged dodatków spoza strict-core.
- [ ] Wejścia/założenia są jawnie wypisane i skończone.

**Pass criterion:** jawny theorem statement + checker summary ze statusem `DISCHARGED`.

### O2. Selector-source identifiability (strict)
- [ ] Kandydat(y) źródła mają jednoznaczny interfejs operacyjny w strict lane.
- [ ] Każdy kandydat daje ten sam kierunek po translacji do wspólnej reprezentacji.
- [ ] Brak niejawnych gauge-fixing hacków.

**Pass criterion:** identifiability status `DISTINGUISHABLE` + brak ukrytych stopni swobody.

### O3. Nonuniqueness exclusion
- [ ] Dla alternatywnych kandydatów wykazano równoważność lub wykluczenie.
- [ ] `nonuniqueness_residual` przechodzi z `OPEN` do `CLOSED`.
- [ ] Counterexample sweep nie znajduje dopuszczalnej kontr-gałęzi.

**Pass criterion:** formal exclusion theorem + replay audit = `PASS`.

### O4. Direction selection theorem
- [ ] Zdefiniowano operator wyboru kierunku w strict lane.
- [ ] Udowodniono niezmienniczość wyboru względem dozwolonych transformacji.
- [ ] Wybór jest stabilny przy perturbacjach w granicach strict tolerance.

**Pass criterion:** theorem status `DISCHARGED`; stability bounds jawne.

### O5. QW-2191 obstruction discharge
- [ ] Główny obstruction bound ma formalne domknięcie (bez luk R5–R11 typu unresolved).
- [ ] Zależności logiczne są acykliczne.
- [ ] Dowód zamyka pełne neighborhood przeszkody, nie tylko fragment.

**Pass criterion:** obstruction status `DISCHARGED_FULL_NEIGHBORHOOD`.

### O6. Independent replay & adversarial audit
- [ ] Niezależne odtworzenie przez drugi przebieg (inna kolejność seed/control where possible).
- [ ] Adversarial checklist: próby rozbicia selektora i wybory kierunku.
- [ ] Wynik zgodny na wszystkich ścieżkach replikacji.

**Pass criterion:** replay `PASS`, divergence count `0`.

### O7. Governance gate: strict claim admissibility
- [ ] Claim label = `STRICT_CORE_CLOSURE` (nie `NON_STRICT`).
- [ ] Brak sprzeczności z guardrails (`QW-2191`, selector guardrail, kernel-split).
- [ ] Export packet rozdziela strict-only od NB-only.

**Pass criterion:** governance gate = `ALLOW_STRICT_QW2191_CLOSURE_CLAIM`.

## Decyzja końcowa
`QW-2191 = CLOSED (strict-only)` jest dopuszczalne **wyłącznie jeśli O1–O7 = PASS**.
W przeciwnym razie status zostaje `NOT_CLOSED`.

## Co z NB?
NB może pozostać domknięte w swoim zakresie i to nie implikuje strict closure.
NB jest torem zgodności/nie-transferu; strict closure wymaga osobnego dowodu wyboru selektora i kierunku.

## Rekomendowany następny uczciwy krok (profesorski)
Uruchomić najpierw **O3 Nonuniqueness exclusion** jako bramkę krytyczną:
1. sformalizować twierdzenie wykluczenia niejednoznaczności,
2. wykonać counterexample sweep,
3. dopiero po `nonuniqueness_residual = CLOSED` przejść do pełnego claimu closure.

## Dla laika (krótko)
Masz kilka sensownych „kompasów” pokazujących kierunek. Ale żeby ogłosić, że teoria jest domknięta, trzeba dowieść, że **wszystkie kompas-y wskazują ten sam kierunek i żaden legalny kompas nie wskazuje innego**. Dopiero wtedy QW-2191 można uczciwie uznać za zamknięte w wersji strict.
