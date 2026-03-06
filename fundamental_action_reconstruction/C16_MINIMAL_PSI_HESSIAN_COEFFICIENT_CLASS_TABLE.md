# C16 Minimal Psi-Hessian Coefficient-Class Table

Status: `C16_EXECUTED_MINIMAL_COEFFICIENT_CLASS_TABLE_NO_FALSE_PASS`
As of: `2026-03-06`

## Cel

Po `C15` najwezszy aktywny blocker brzmial:

- `C15_B1 := no_explicit_coefficient_filled_canonical_Psi_x_Psi_block_H_PsiPsi_for_evaluating_the_control_pullback`

`C16` nie probuje udawac, ze strict core eksportuje juz pelna,
wyczerpujaca macierz wspolczynnikow `12 x 12` dla canonical `H_PsiPsi`.

`C16` robi cos wezszejszego:
- sprawdza, czy strict core zawiera juz minimalny coefficient-class table
  dla reprezentatywnych wierszy `Psi` sektora,
- i czy blocker da sie zawezic do braku exhaustive canonical coefficient table
  oraz braku restriction do orientation slice.

## Polityka zrodel

### Strict-admissible support

1. `QW-2164`
   - `sample_linearized_eom_eta0`,
   - `n_psi_fields = 12`.
2. `QW-2166`
   - `sample_eom_eta0`,
   - `sample_eom_eta6`,
   - `n_psi_fields = 12`,
   - `linear_operator_matrix_matches_canonical_hessian=True`.
3. `QW-2180`
   - exact operator/Hessian identification.
4. `C12`
   - coefficient classes already isolated at packet level.
5. `C15`
   - control-only pullback packet `M_control = T_control^T H_PsiPsi T_control`.
6. `A10`
   - anti-overclaim boundary.

## Pytanie audytowe

Czy wolno juz powiedziec, ze:
- strict core ma juz przynajmniej reprezentatywny coefficient-class table
  dla canonical `H_PsiPsi`,
- bo sample rows `eta0` i `eta6` jawnie eksportuja klasy wpisow
  diagonalnych i off-diagonalnych,
- a wiec blocker nie dotyczy juz braku jakiegokolwiek coefficient filling,
  tylko braku exhaustive canonical table dla calego `12 x 12` bloku?

## Co daje `QW-2164`

`QW-2164` daje dla `eta0`:
- off-diagonal kernel-mixing entries typu
  - `(K_0_j + K_j_0)/2` w warstwie sample/proxy-expanded,
- diagonalny pakiet lokalny:
  - `3*g4_psi0*vpsi0**2`,
  - `5*g6_psi0*vpsi0**4`,
  - `2*gY0*vphi**2`,
  - `m2_psi0`,
- kinetic identity term:
  - `Derivative(eta0(x), (x, 2))`.

To jest juz coefficient-class row dla reprezentatywnego `Psi` indeksu.

## Co daje `QW-2166`

`QW-2166` wzmacnia to do exhaustive layer:
- `eta0` ma jawny row z symetryzowanym mixingiem
  `(K_0_j + K_j_0)/2`,
- `eta6` daje drugi jawny row tego samego typu,
- diagonalny schemat dla `eta6` ma identyczna klase:
  - `3*g4_psi6*vpsi6**2`,
  - `5*g6_psi6*vpsi6**4`,
  - `2*gY6*vphi**2`,
  - `m2_psi6`,
- operator matrix matches canonical Hessian.

To znaczy:
- coefficient-class filling nie jest juz pustym placeholderem,
- strict core ma juz dwa reprezentatywne rows canonical `Psi` blocku.

## Minimalna tabela po `C16`

Najmocniejsza uczciwa tabela klas wspolczynnikow brzmi:

| Row seed | Off-diagonal class | Diagonal class | Kinetic class |
|---|---|---|---|
| `eta0` | `(K_0_j + K_j_0)/2` | `3*g4_psi0*vpsi0**2 + 5*g6_psi0*vpsi0**4 + 2*gY0*vphi**2 + m2_psi0` | `d^2/dx^2` |
| `eta6` | `(K_6_j + K_j_6)/2` | `3*g4_psi6*vpsi6**2 + 5*g6_psi6*vpsi6**4 + 2*gY6*vphi**2 + m2_psi6` | `d^2/dx^2` |

To nie jest jeszcze pelna tabela dla wszystkich `12` indeksow `Psi`.

## Wynik `C16`

`C16` ustala:
- `C15_B1` nie dotyczy juz braku jakiegokolwiek coefficient-class filling,
- strict core ma juz minimalny coefficient-class table
  dla reprezentatywnych rows canonical `H_PsiPsi`,
- control-only pullback z `C15` ma juz packet-ready input classes,
  nawet jesli nie ma jeszcze exhaustive canonical table.

## Redukcja frontu po `C16`

Po `C15` mielismy:

- `C15_B1 := no_explicit_coefficient_filled_canonical_Psi_x_Psi_block_H_PsiPsi_for_evaluating_the_control_pullback`
- `C15_B2 := no_explicit_restriction_from_M_control_to_the_candidate_orientation_slice`

Po `C16` najuczciwiej zapisac to weziej jako:

- `C16_B1 := no_exhaustive_index_complete_coefficient_table_for_the_canonical_12x12_Psi_x_Psi_block_H_PsiPsi`
- `C16_B2 := no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice`

To jest realny postep redukcyjny:
- minimalna tabela klas wspolczynnikow juz istnieje,
- otwarty brak dotyczy tabeli exhaustive i restriction do slice.

## Macierz wyniku

| Pytanie | Status po C16 | Uwagi |
|---|---|---|
| representative coefficient-class rows exist | `present_partial` | `eta0`, `eta6` |
| coefficient classes cover diagonal/off-diagonal/kinetic structure | `present_partial` | nadal bez exhaustive full matrix |
| exhaustive `12 x 12` coefficient table exists | `not_shown` | nadal brak |
| restriction to candidate orientation slice exists | `not_shown` | nadal brak |
| discharge of `C15_B1` | `reduced_not_closed` | blocker tylko zawężony |

## Czego `C16` nie ustala

`C16` nie ustala:
- ze cala canonical macierz `H_PsiPsi` zostala wyeksportowana,
- ze kazdy indeks `0..11` ma juz jawny row w coefficient table,
- ze `M_control` zostal juz policzony liczbowo lub symbolicznie do konca,
- ze orientation-slice restriction istnieje,
- ze `C15_B1` ma PASS,
- ze `C15_B2` ma PASS.

## Anti-overclaim

`C16` nie twierdzi, ze:
- dwa sample rows wystarczaja do pelnego theorem-level block matching,
- reprezentatywna tabela klas jest rownowazna exhaustive canonical table,
- restriction do orientation slice jest blisko domkniecia.

## Produkt etapu

- szesnasty krok trzeciego mikrocyklu,
- minimalna coefficient-class table dla canonical `H_PsiPsi`,
- zawężenie `C15_B1` do braku exhaustive index-complete table,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `C17`:
- sprawdzic, czy strict core pozwala juz na finite index-complete stencil
  dla wszystkich `12` rows canonical `H_PsiPsi`,
- albo jawnie potwierdzic, ze exhaustive coefficient export nadal nie jest obecny.
