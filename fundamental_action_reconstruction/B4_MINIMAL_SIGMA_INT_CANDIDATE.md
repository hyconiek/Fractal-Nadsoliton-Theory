# B4 Minimal Sigma_int Candidate

Status: `B4_EXECUTED_HYBRID_SIGMA_INT_CANDIDATE_IDENTIFIED_STRICT_DISCHARGE_PENDING`
As of: `2026-03-06`

## Cel

`B4` podejmuje tylko `B3_O1`:
- zdefiniowac minimalny kandydat `sigma_int`,
- bez twierdzenia, ze selector zostal juz wyprowadzony theorem-level.

## Polityka zrodel

### Strict-admissible support

1. `QW-2206`
   - local topological protection layer integrated.
2. `B3`
   - packet `B3_O1..B3_O5`.
3. `A10`
   - anti-overclaim boundary.

### Hybrid support

1. `QW-1622`
   - FR quantization gives binary sign structure on the nontrivial loop class.

## Kandydat

Niech:
- `gamma_pi1` oznacza generator nietrywialnej klasy petli w `pi_1(C)=Z_2`
  dla lokalnego sektora topologicznego jednego nadsolitonu,
- `chi_FR(gamma_pi1)` oznacza znak FR przypisany tej klasie.

Wtedy minimalny kandydat definiujemy jako:

`sigma_int_candidate := chi_FR(gamma_pi1) in {+1, -1}`.

W lokalnym sektorze `B~1` z warstwa FR:
- kandydat przyjmuje naturalnie wartosc `-1`.

## Dlaczego to jest dobry kandydat

1. Jest binarny
   - pasuje do najmniejszej mozliwej dodatkowej struktury.
2. Jest wewnetrzny
   - nie pochodzi z zewnetrznego convention token.
3. Jest topologiczny
   - nie jest zwyklym fit-parametrem.
4. Jest zgodny z ontologia jednego nadsolitonu
   - bo opiera sie na strukturze samego obiektu i jego przestrzeni konfiguracji.

## Co `B4` rzeczywiscie ustala

`B4` ustala:
- najprostszy sensowny kandydat `sigma_int` nie jest juz abstrakcyjnym placeholderem,
- tylko konkretnym obiektem:
  - FR-sign of the nontrivial loop class.

To jest realny postep, bo `B3_O1` ma teraz jeden kanoniczny kandydat zamiast wielu nieporownywalnych hipotez.

## Czego `B4` jeszcze nie ustala

`B4` nie ustala:
- ze `sigma_int_candidate` jest juz strict-derived,
- ze `sigma_int_candidate` jest stabilny na calej potrzebnej rodzinie deformacji,
- ze jest gauge-safe w sensie potrzebnym dla `QW-2191`,
- ze istnieje juz most:
  - `sigma_int_candidate -> theta selector`,
- ze axiom-free uniqueness zostala zamknieta.

## Status obligacji `B3_O1..B3_O5` po B4

| Obligacja | Status po B4 | Uwagi |
|---|---|---|
| `B3_O1` define internal datum | `candidate identified` | FR-sign candidate |
| `B3_O2` deformation/gauge stability | `open` | brak dowodu |
| `B3_O3` map datum to selector | `open` | brak mostu do `theta` |
| `B3_O4` compatibility with mode scaffold | `open` | brak theorem package |
| `B3_O5` anti-overclaim closure test | `pending after O2..O4` | jeszcze za wczesnie |

## Anti-overclaim

`B4` nie twierdzi, ze:
- FR sign automatycznie rozladowuje `QW-2191`,
- `sigma_int_candidate` jest juz the one true selector,
- lokalna warstwa topologiczna automatycznie daje globalny uniqueness theorem.

`B4` twierdzi tylko:
- istnieje jeden minimalny, kanoniczny kandydat `sigma_int`.

## Produkt etapu

- czwarty krok drugiego cyklu,
- kanoniczny kandydat dla `B3_O1`,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `B5`:
- sprawdzic `B3_O2`:
  - czy `sigma_int_candidate` jest stabilny na dopuszczalnych deformacjach
    i nie znika po poprawnym quotient przez gauge / parametryzacje.
