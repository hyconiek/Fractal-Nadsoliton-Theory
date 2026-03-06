# C13 Mode-Basis Control Index-Set Audit

Status: `C13_EXECUTED_MODE_BASIS_CONTROL_INDEX_SET_REDUCTION_NO_FALSE_PASS`
As of: `2026-03-06`

## Cel

Po `C12` najwezszy blocker brzmial:

- `C12_B1 := no_explicit_assembled_Psi_x_Psi_submatrix_and_no_coefficient_table_for_a_chosen_index_set_ready_for_matching_against_qw2186`

`C13` nie probuje udawac, ze strict core wybral juz kanoniczny `Psi`-index-set w carrierze Hessianu.

`C13` robi cos wezszejszego:
- sprawdza, czy strict core zawiera juz przynajmniej kontrolny index-set w bazie modowej,
- i czy blocker da sie zawezic do braku transportu z tej bazy modowej do canonical `Psi` basis.

## Polityka zrodel

### Strict-admissible support

1. `QW-2190`
   - deterministic real Fourier mode basis,
   - `SU(3): [e0,c1,s1]`,
   - `SU(2): [c2,s2]`,
   - full physical uniqueness of mode-index assignment remains open.
2. `C3`
   - reference pair candidates `(c1,s1)` and `(c2,s2)`.
3. `C7`
   - class-level schema `mode pair -> orientation slice`.
4. `C12`
   - minimal `Psi-sector block` extraction packet.
5. `A10`
   - anti-overclaim boundary.

## Pytanie audytowe

Czy wolno juz powiedziec, ze:
- strict core ma przynajmniej deterministyczne kontrolne zbiory indeksow w bazie modowej,
- a wiec blocker nie dotyczy juz wyboru index-setu w ogole,
- tylko transportu:
  - `mode-basis control index-set -> canonical Psi index-set for Hessian extraction`?

## Co daje `QW-2190`

`QW-2190` ustala deterministycznie:
- `pair1 := (c1,s1)`,
- `pair2 := (c2,s2)`,
- bez skanu i bez retune,
- ale z jawnie otwarta pelna fizyczna unikalnoscia przypisania indeksow modow.

To daje dwa uczciwe control index-sets:

- `I_mode_1 := {c1, s1}`
- `I_mode_2 := {c2, s2}`

Nie sa to jeszcze canonical `Psi` index-sets.

## Wynik `C13`

`C13` ustala:
- strict core ma juz jawne control index-sets w bazie modowej,
- brak z `C12` nie dotyczy juz wyboru `I` w ogole,
- lecz braku jawnego transportu od `I_mode_1` / `I_mode_2` do canonical `Psi` basis carrieru Hessianu.

## Redukcja frontu po `C13`

Po `C12` mielismy:

- `C12_B1 := no_explicit_assembled_Psi_x_Psi_submatrix_and_no_coefficient_table_for_a_chosen_index_set_ready_for_matching_against_qw2186`

Po `C13` najuczciwiej zawezic to do:

- `C13_B1 := no_explicit_transport_from_the_deterministic_mode_basis_control_index_set_to_a_canonical_Psi_index_set_inside_the_exhaustive_Hessian_carrier`
- `C13_B2 := no_assembled_Psi_x_Psi_submatrix_after_such_transport`

To jest realny postep redukcyjny:
- index-set selection exists in control mode basis,
- otwarty pozostaje transport i assembled submatrix.

## Macierz wyniku

| Pytanie | Status po C13 | Uwagi |
|---|---|---|
| deterministic control index-set in mode basis exists | `present` | `I_mode_1`, `I_mode_2` |
| canonical Psi index-set chosen | `not_shown` | nadal brak |
| transport mode basis -> Psi basis | `not_shown` | nadal brak |
| assembled Psi x Psi submatrix after transport | `not_shown` | nadal brak |
| discharge of `C12_B1` | `split_not_closed` | blocker tylko rozbity |

## Czego `C13` nie ustala

`C13` nie ustala:
- ze `I_mode_1` albo `I_mode_2` jest juz physical canonical choice,
- ze transport do `Psi` basis zostal juz skonstruowany,
- ze assembled submatrix istnieje,
- ze `C12_B1` ma PASS,
- ze `C9_B2` ma PASS.

## Anti-overclaim

`C13` nie twierdzi, ze:
- control index-set w bazie modowej rozwiazuje problem block matchingu,
- `QW-2190` juz daje kanoniczne przypisanie do carrieru Hessianu,
- theorem-level uniqueness closure jest blisko.

## Produkt etapu

- trzynasty krok trzeciego mikrocyklu,
- jawny control index-set w bazie modowej,
- rozbicie `C12_B1` na brak transportu i brak assembled submatrix,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `C14`:
- sprobowac minimalnego transport schema `I_mode -> Psi basis`,
- albo jawnie potwierdzic, ze strict core nadal nie eksportuje takiego mostu.
