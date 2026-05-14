# P1662 / S612 Strict Kernel -> Full Lagrangian -> EOM Bidirectional Witness Packet (No Legacy Bridge)

Status: `P1662_EXECUTED_STRICT_FULL_CHAIN_WITNESS_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel

Wykonać następny uczciwy krok strict-only w torze:

`K_strict -> współczynniki -> pełny L_SM+L_GR -> równania ruchu`

i jawnie pokazać kierunek odwrotny (`EOM constraints -> coefficient consistency`) bez fałszywego pass.

## Zakres fizyczny

W pakiecie eksportujemy *operacyjny pełny ansatz* gęstości lagranżianu (SM+GR) w ścisłej postaci roboczej:

- sektor grawitacyjny: `sqrt(-g) * (Mpl2/2 * R + cR2 * R^2 + cRic2 * R_{mu nu}R^{mu nu})`,
- sektor cechowania: `-1/4 sum_i Z_i F_i^2`,
- sektor Higgsa: `|D H|^2 - muH2 |H|^2 - lambdaH |H|^4`,
- sektor fermionowy: `i psi_bar gamma^mu D_mu psi - y psi_bar H psi + h.c.`,
- sektor mieszania strict: `xiHR |H|^2 R + chiRG R tr(F^2)`.

Bez użycia bridge do legacy.

## Warstwa kernel -> współczynniki

Stosujemy jawny strict-only map-maker (operacyjny, nie ontologiczny):

- wejście kernela: `(omega, phi, beta, eta, A)`,
- wyjście współczynników EFT: `(Mpl2, cR2, cRic2, Z_i, muH2, lambdaH, y, xiHR, chiRG)`,
- wymóg: dodatniość wybranych norm i brak sprzeczności znaków energetycznych.

## Warstwa Lagrangian -> EOM

Eksport CAS (`p1662_strict_full_chain_symbolic_witness.py`) wykonuje:

1. budowę symboliczną `L_total`,
2. równanie Eulera-Lagrange'a dla pola skalarnego zastępczego (`H -> h`),
3. równanie pola dla U(1)-proxy,
4. testy spójności znaku i redukcji do limitu swobodnego.

## Kierunek odwrotny (EOM -> współczynniki)

Wymuszamy kontrolę odwrotną:

- z EOM skalarnego odzyskujemy kombinacje `muH2, lambdaH, xiHR`,
- z EOM gauge odzyskujemy `Z_A` i warunek źródłowy,
- status PASS tylko gdy residuale symboliczne redukują się do zera.

## QG: rygory ToE (bez false pass)

Pakiet **nie** ogłasza zamknięcia QG; tylko ustanawia macierz zobowiązań strict-core:

- renormalizacja: wymagany eksport liczenia kontrterminów i beta-funkcji,
- unitarność: wymagane widmo bez duchów dla operatorów wyższej krzywizny,
- background independence: wymagany wariacyjny zapis i test niezależności od tła.

Brak któregoś eksportu => `OPEN_OBLIGATION`, nie PASS.

## Rekomendowany następny uczciwy krok

Uruchomić `S613`: jawny moduł kontrterminów 1-loop dla `R^2/R_{mu nu}^2` + test unitarności propagatora w sektorze spin-2.

## Omówienie dla laika

To jest jak budowa silnika z planu technicznego:

1. kernel strict to surowe ustawienia,
2. współczynniki to konkretne części silnika,
3. lagranżian to pełna instrukcja działania,
4. równania ruchu to test, czy silnik faktycznie pracuje.

W tym kroku pokazujemy, że instrukcja i test pasują do siebie w obie strony,
ale uczciwie zaznaczamy, że pełna „certyfikacja bezpieczeństwa” dla grawitacji kwantowej
(renormalizacja + unitarność + niezależność od tła) jest jeszcze otwarta.
