# C28 Local Orbit-Frame Quotient Schema

Status: `C28_EXECUTED_LOCAL_ORBIT_FRAME_QUOTIENT_SCHEMA_NO_FALSE_PASS`
As of: `2026-03-06`

## Cel

Po `C27` najwezszy aktywny blocker brzmi:

- `C27_B1 := no_explicit_control_coordinate_realization_of_the_zero_mode_quotient_candidate_on_the_control_pullback_orbit_family`

`C28` nie probuje udawac, ze strict core ma juz globalny, jawnie
wyeksportowany quotient operator na calej rodzinie `M_control`.

`C28` robi cos wezszejszego:
- sprawdza, czy w lokalnych control coordinates
  `span(c1,s1)` albo `span(c2,s2)`
  istnieje juz packet-ready schema lokalnej dekompozycji:
  - kierunek orbit tangent,
  - kierunek poprzeczny / mismatch.

## Polityka zrodel

### Strict-admissible support

1. `C4`
   - lokalna geometria orbity `O(2)`,
   - `u(theta)`, `v(theta)`,
   - `Delta u`, `Delta v`,
   - kinematyczny split na orbicie.
2. `C14`
   - control transport schema,
   - control basis `(c1,s1,c2,s2)`.
3. `C15`
   - control-only pullback carrier `M_control`.
4. `C27`
   - quotient candidate class po odjeciu modow zerowych.
5. `A10`
   - anti-overclaim boundary.

## Pytanie audytowe

Czy wolno juz powiedziec, ze strict core ma packet-ready
**lokalny** control-coordinate quotient schema w nastepujacym sensie:

- dla kazdej pary `(c_i,s_i)` istnieje lokalna rama orbity,
- jedna os jest tangent do orbity,
- druga os jest poprzeczna i naturalnie odpowiada lokalnemu mismatch
  direction,
- wiec brak jawnego quotient map nie jest juz brakiem jakiejkolwiek
  control-coordinate realizacji, tylko brakiem jej jawnej serializacji
  i sklejenia globalnego?

## Co daje `C4`

`C4` ustala lokalny obraz orbity:

- `u(theta) = cos(theta) c_ref + sin(theta) s_ref`,
- `v(theta) = -sin(theta) c_ref + cos(theta) s_ref`.

Na tej orbicie:
- tangent direction jest dana przez pochodna po `theta`,
- mismatch direction jest opisana przez `Delta u(theta)` i `Delta v(theta)`,
- cross-term zanika identycznie.

To nie daje jeszcze globalnego projektora.
Ale daje lokalny orbit-frame schema.

## Lokalna control-coordinate realizacja

W `span(c_ref,s_ref)` mozna zapisac lokalnie:

- tangent axis:
  - `tau(theta) := -sin(theta) c_ref + cos(theta) s_ref`,
- mismatch / radial axis:
  - `nu(theta) := cos(theta) c_ref + sin(theta) s_ref - c_ref`,
    albo rownowazna znormalizowana wersja poprzeczna.

Najuczciwszy wniosek:
- strict core ma juz lokalny control-coordinate kandydat
  na quotient split:
  - quotient by orbit tangent = retain transverse mismatch direction.

To jest schema lokalne, zalezne od punktu orbity `theta`.

## Co daje `C14` i `C15`

`C14` daje:
- control basis osadzona w canonical carrierze `Psi`.

`C15` daje:
- control-only pullback carrier `M_control`.

Razem oznacza to:
- lokalny orbit-frame schema nie jest juz odklejony od aktualnego nośnika,
- ale nadal nie jest wyeksportowany jako jawny operator dzialajacy
  na calej rodzinie `M_control`.

## Najmocniejszy uczciwy wniosek po `C28`

Po zlozeniu `C4 + C14 + C15 + C27` najuczciwiej zapisac:

- istnieje packet-ready **lokalny** control-coordinate quotient schema,
- jego tresc to:
  - identify local orbit tangent direction,
  - quotient it out,
  - keep the transverse mismatch direction as reduced control-plane candidate.

To jest realny postep.

Nie jest to jeszcze:
- globalny, jawny, serialized projector,
- ani finalna extraction map do dwuwymiarowej orientation slice.

## Redukcja frontu po `C28`

Po `C27` mielismy:

- `C27_B1 := no_explicit_control_coordinate_realization_of_the_zero_mode_quotient_candidate_on_the_control_pullback_orbit_family`
- `C26_B2 := no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane`

Po `C28` najuczciwiej zapisac pierwszy blocker weziej jako:

- `C28_B1 := no_explicit_serialized_local_orbit_frame_projector_formula_or_global_gluing_rule_for_the_control_coordinate_quotient_candidate`

Drugi blocker pozostaje bez zmian:

- `C26_B2 := no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane`

To jest realny postep redukcyjny:
- local control-coordinate realization candidate juz jest,
- brak dotyczy serializacji / global gluing.

## Macierz wyniku

| Pytanie | Status po C28 | Uwagi |
|---|---|---|
| local orbit-frame schema in `span(c_i,s_i)` exists | `present_partial` | `C4` |
| local control-coordinate quotient candidate exists | `present_partial` | tangent vs transverse |
| explicit serialized projector formula exists | `not_shown` | nadal brak |
| global gluing rule exists | `not_shown` | nadal brak |
| final basis-level slice extraction exists | `not_shown` | nadal brak |
| discharge of `C27_B1` | `reduced_not_closed` | tylko zawężenie |

## Czego `C28` nie ustala

`C28` nie ustala:
- ze jawny projector zostal zapisany jako macierz lub operator,
- ze quotient schema jest globalnie sklejony po calej rodzinie orbit,
- ze local orbit-frame jest juz canonical theorem-level,
- ze finalna orientation slice zostala wydobyta,
- ze `C27_B1` ma PASS.

## Anti-overclaim

`C28` nie twierdzi, ze:
- local orbit-frame schema jest juz gotowa restriction map,
- quotient problem jest rozwiazany globalnie,
- orientation-slice closure jest blisko theorem-level PASS.

## Produkt etapu

- dwudziesty ósmy krok trzeciego mikrocyklu,
- local control-coordinate quotient schema,
- zawężenie `C27_B1` do braku serialized projector / global gluing,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `C29`:
- sprawdzic, czy strict core ma juz packet-ready serialized formula
  dla lokalnego orbit-frame projector,
- albo jawnie potwierdzic, ze pozostaje tylko ten brak i finalny slice extraction.
