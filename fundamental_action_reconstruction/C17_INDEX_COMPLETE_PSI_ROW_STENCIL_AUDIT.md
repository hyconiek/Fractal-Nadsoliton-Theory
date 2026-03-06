# C17 Index-Complete Psi Row Stencil Audit

Status: `C17_EXECUTED_INDEX_COMPLETE_ROW_STENCIL_NO_FALSE_PASS`
As of: `2026-03-06`

## Cel

Po `C16` najwezszy aktywny blocker brzmial:

- `C16_B1 := no_exhaustive_index_complete_coefficient_table_for_the_canonical_12x12_Psi_x_Psi_block_H_PsiPsi`

`C17` nie probuje udawac, ze strict core eksportuje juz pelna,
row-by-row tablice wspolczynnikow dla calego canonical bloku `H_PsiPsi`.

`C17` robi cos wezszejszego:
- sprawdza, czy strict core zawiera juz index-complete **row stencil schema**
  dla wszystkich `12` wierszy `Psi`,
- i czy blocker da sie zawezic do braku jawnej instancjacji
  wspolczynnikow dla kazdego `i=0..11`,
  a nie do braku samego wzoru wiersza.

## Polityka zrodel

### Strict-admissible support

1. `QW-2163`
   - canonical action `12xPsi + Phi`,
   - explicit index mixing `K_{i,j}`.
2. `QW-2165`
   - exhaustive canonical EoM for all `13` fields,
   - locality/self/Yukawa/kernel-mixing verified exhaustively.
3. `QW-2166`
   - exhaustive canonical Hessian for all `13` fields,
   - `n_psi_fields = 12`,
   - `linear_operator_matrix_matches_canonical_hessian=True`.
4. `QW-2180`
   - exact operator/Hessian identification.
5. `C16`
   - representative coefficient-class rows for `eta0` and `eta6`.
6. `A10`
   - anti-overclaim boundary.

## Pytanie audytowe

Czy wolno juz powiedziec, ze:
- strict core ma juz skonczony, index-complete stencil wiersza dla kazdego
  `Psi_i`,
- bo exhaustive canonical EoM/Hessian obejmuja wszystkie `12` pola `Psi`,
- a wiec blocker nie dotyczy juz samej **postaci** wiersza,
  tylko braku jawnego coefficient export dla wszystkich `i=0..11`?

## Co daje `QW-2163`

`QW-2163` ustala canonical action:
- `12xPsi + Phi`,
- explicit index mixing `K_{i,j}`,
- lokalny carrier dzialania dla wszystkich pol.

To daje wspolny action-level szablon sprzezen miedzy `Psi_i` i `Psi_j`.

## Co daje `QW-2165`

`QW-2165` wzmacnia ten obraz do poziomu exhaustive:
- EoM sa policzone dla wszystkich `12` pol `Psi` oraz `Phi`,
- locality, self-interaction, Yukawa i kernel-mixing sa sprawdzone
  exhaustively, a nie tylko na jednym sample row.

To znaczy:
- wzor wiersza `Psi_i` nie jest juz ograniczony do `eta0` lub `eta6`,
- istnieje finite row schema dla wszystkich `i = 0..11`.

## Co daje `QW-2166`

`QW-2166` przenosi to na canonical Hessian:
- `n_psi_fields = 12`,
- exhaustive canonical Hessian dla wszystkich `13` pol,
- exact matrix/operator consistency.

Najslabszy uczciwy stencil dla kazdego `Psi_i` brzmi:

- diagonal:
  - `3*g4_psii*vpsii**2 + 5*g6_psii*vpsii**4 + 2*gYi*vphi**2 + m2_psii`
- off-diagonal `Psi-Psi`:
  - `(K_i_j + K_j_i)/2` dla `j != i`
- kinetic identity:
  - `d^2/dx^2`

Pelny row moze miec tez sprzezenie do `Phi`,
ale ono nie nalezy do canonical bloku `H_PsiPsi`.

## Wynik `C17`

`C17` ustala:
- strict core ma juz index-complete row stencil schema
  dla wszystkich `12` wierszy canonical `H_PsiPsi`,
- `C16_B1` nie dotyczy juz braku wzoru wiersza,
- tylko braku jawnej row-by-row instancjacji wspolczynnikow
  dla wszystkich `i=0..11`.

## Redukcja frontu po `C17`

Po `C16` mielismy:

- `C16_B1 := no_exhaustive_index_complete_coefficient_table_for_the_canonical_12x12_Psi_x_Psi_block_H_PsiPsi`
- `C16_B2 := no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice`

Po `C17` najuczciwiej zapisac to weziej jako:

- `C17_B1 := no_explicit_row_by_row_export_instantiating_the_index_complete_Psi_row_stencil_for_all_i_0_to_11`
- `C17_B2 := no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice`

To jest realny postep redukcyjny:
- index-complete row stencil jest juz obecny,
- otwarty brak dotyczy explicit row-by-row exportu i restriction do slice.

## Macierz wyniku

| Pytanie | Status po C17 | Uwagi |
|---|---|---|
| index-complete row stencil schema exists | `present_partial` | dla wszystkich `Psi_i`, `i=0..11` |
| exhaustive row-by-row coefficient export exists | `not_shown` | nadal brak |
| exhaustive `12 x 12` canonical table exists | `not_shown` | nadal brak |
| restriction to candidate orientation slice exists | `not_shown` | nadal brak |
| discharge of `C16_B1` | `reduced_not_closed` | blocker tylko zawężony |

## Czego `C17` nie ustala

`C17` nie ustala:
- ze kazdy row `i=0..11` jest juz osobno wyeksportowany,
- ze cala canonical macierz `H_PsiPsi` jest juz zapisana wspolczynnik po wspolczynniku,
- ze control pullback jest juz policzony,
- ze orientation-slice restriction istnieje,
- ze `C16_B1` ma PASS,
- ze `C16_B2` ma PASS.

## Anti-overclaim

`C17` nie twierdzi, ze:
- exhaustive EoM/Hessian automation jest rownowazna pelnemu coefficient exportowi,
- index-complete stencil oznacza juz gotowa macierz `12 x 12`,
- restriction do candidate orientation slice jest blisko domkniecia theorem-level.

## Produkt etapu

- siedemnasty krok trzeciego mikrocyklu,
- index-complete row stencil schema dla wszystkich `12` pol `Psi`,
- zawężenie `C16_B1` do braku jawnego row-by-row exportu,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `C18`:
- sprawdzic, czy strict core pozwala juz na finite row-by-row export packet
  dla wszystkich `i=0..11`,
- albo jawnie potwierdzic, ze exhaustive coefficient export nadal nie jest obecny.
