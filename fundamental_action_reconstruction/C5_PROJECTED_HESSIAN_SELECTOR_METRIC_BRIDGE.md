# C5 Projected Hessian Selector-Metric Bridge

Status: `C5_EXECUTED_PROJECTED_HESSIAN_BRIDGE_CONDITIONAL_LOCAL_IDENTIFICATION_PENDING_NO_FALSE_PASS`
As of: `2026-03-06`

## Cel

Po `C4` wiemy juz, ze jesli istnieje lokalna dodatnia metryka mismatch na kandydackiej plaszczyznie orientacji, to forma `J_ab(theta)` redukuje sie do:

`J_ab(theta)=2(a+b)(1-cos theta)`.

`C5` sprawdza cos mocniejszego:
- czy ten lokalny koszt da sie uczciwie utozsamic nie z dowolna metryka "z zewnatrz",
- tylko z naturalnym kandydatem pochodzacym z drugiej wariacji / projected Hessian,
- i ile z tego wynika juz teraz bez falszywego PASS.

## Polityka zrodel

### Strict-admissible support

1. `A3`
   - operator drugiej wariacji,
   - projection-before-claim discipline.
2. `A7`
   - branch-scope positivity package.
3. `QW-2190`
   - deterministic mode scaffold i kandydackie pary degenerowane.
4. `QW-2191`
   - obstruction theorem dla fizycznej unikalnosci.
5. `C3`
   - reference-pair candidate.
6. `C4`
   - kinematyczna redukcja kosztu mismatch.
7. `A10`
   - anti-overclaim boundary.

### Jawnie poza rdzeniem C5

`C5` nie zaklada:
- ze projected Hessian na kandydackiej plaszczyznie orientacji zostal juz explicite wyliczony,
- ze jego dodatniosc na tej plaszczyznie zostala juz certyfikowana,
- ze orientacyjna plaszczyzna jest juz podniesiona do physical selector datum.

## Pytanie audytowe

Czy wolno powiedziec, ze jesli:
1. istnieje projekcja drugiej wariacji na kandydacka plaszczyzne orientacji,
2. ta projekcja daje lokalna symetryczna forme kwadratowa,
3. forma ta ma dodatni certyfikat in-scope,

to selector family `J_ab` jest juz naturalna forma orbitalna projected Hessianu, a nie arbitralnym kosztem dodanym recznie?

## Ustawienie

Po `C3` bierzemy kandydacka pare:
- `(c_ref, s_ref) := (c1, s1)` albo `(c2, s2)`.

Na orbicie `O(2)`:
- `u(theta) = cos(theta) c_ref + sin(theta) s_ref`
- `v(theta) = -sin(theta) c_ref + cos(theta) s_ref`
- `Delta u = u(theta) - c_ref`
- `Delta v = v(theta) - s_ref`

## Naturalny kandydat pochodzacy z drugiej wariacji

Jesli projected second variation istnieje na tej kandydackiej plaszczyznie, to najogolniejszy lokalny symetryczny koszt kwadratowy ma postac:

`Q_H(theta) = a <Delta u,Delta u> + 2c <Delta u,Delta v> + b <Delta v,Delta v>`

z:
- `a,b,c` lokalnymi wspolczynnikami projected Hessianu,
- dodatniosc w strict scope oznacza, ze ta forma ma certyfikat dodatni na rozwazanej plaszczyznie.

## Wynik kinematyczny

Na orbicie `O(2)` z `C4` mamy identycznie:

- `<Delta u,Delta u> = 2(1-cos theta)`
- `<Delta v,Delta v> = 2(1-cos theta)`
- `<Delta u,Delta v> = 0`

W konsekwencji:

`Q_H(theta)=2(a+b)(1-cos theta)`

niezaleznie od wspolczynnika mieszanego `c`.

To jest realne wzmocnienie wzgledem `C4`:
- nie trzeba juz nawet zakladac z gory diagonalnosci,
- wystarcza standardowa lokalna symetryczna forma kwadratowa pochodzaca z projected Hessianu.

## Co `C5` rzeczywiscie ustala

`C5` ustala:
- jesli selector metric ma pochodzic z projected second variation,
  to na orbicie `O(2)` jego orbitalna forma automatycznie redukuje sie do rodziny `J_ab`,
- selector family z `QW-2192/2193` jest zgodna z naturalnym projected-Hessian picture,
- o ile tylko taka projekcja i jej dodatniosc zostana rzeczywiscie zidentyfikowane.

## Czego `C5` nie ustala

`C5` nie ustala:
- ze projected Hessian na kandydackiej plaszczyznie orientacji zostal juz explicite wyciagniety z strict core,
- ze dodatniosc tej projekcji zostala juz udowodniona,
- ze `A3` i `A7` juz razem rozladowuja selector origin,
- ze `C2_B2` ma PASS,
- ze uniqueness jest theorem-level closed.

## Macierz wyniku

| Pytanie | Status po C5 | Uwagi |
|---|---|---|
| standard symmetric projected-Hessian form reduces to `J_ab` on orbit | `derived_conditionally` | cross-term znika identycznie |
| need for diagonal assumption from `C4` | `removed` | wystarcza lokalna symetryczna forma kwadratowa |
| explicit projected Hessian on candidate orientation plane | `not_shown` | strict core tego jeszcze nie eksportuje |
| positivity certificate for that projection | `not_shown` | `A7` nie daje tego na tej plaszczyznie |
| discharge of `C2_B2` | `reduced_not_closed` | brak internal identification pozostaje |

## Redukcja frontu po C5

Po `C4` mielismy:
- `C4_B1 := no_internal_identification_of_the_physical_positive_local_metric_on_candidate_orientation_plane`

Po `C5` najuczciwiej zapisac frontier jako:

- `C5_B1 := no_explicit_projected_second_variation_with_strict_scope_positivity_certificate_on_candidate_orientation_plane`

To jest wezszy blocker, bo:
- wiadomo juz, jaka forma projected Hessianu musi dawac selector family na orbicie,
- otwarte pozostaje juz tylko:
  - realne wyciecie tej projekcji z rdzenia,
  - i realny certyfikat dodatniosci na tej plaszczyznie.

## Anti-overclaim

`C5` nie twierdzi, ze:
- projected Hessian zostal znaleziony,
- `A3` lub `A7` juz rozladowaly selector origin,
- `QW-2192/2193` zostaly zinternalizowane,
- `C2_B2` ma PASS,
- `QW-2191` zostalo rozladowane.

## Produkt etapu

- piaty krok trzeciego mikrocyklu,
- conditional bridge:
  - `projected Hessian -> selector family orbital form`,
- usuniecie zbednego zalozenia diagonalnosci,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `C6`:
- sprawdzic, czy strict core zawiera choc packet-ready kandydat projekcji drugiej wariacji na kandydacka plaszczyzne orientacji,
- albo jawnie potwierdzic, ze tego eksportu nadal nie ma.
