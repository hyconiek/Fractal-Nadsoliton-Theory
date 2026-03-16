# C15 Control-Only Pullback Submatrix Packet

Status: `C15_EXECUTED_CONTROL_ONLY_PULLBACK_PACKET_NO_FALSE_PASS`
As of: `2026-03-16`

## Cel

Po `C14` najwezszy aktywny blocker brzmial:

- `C14_B2 := no_assembled_Psi_x_Psi_submatrix_after_adopting_the_control_transport_schema`

`C15` nie probuje udawac, ze strict core ma juz **fizycznie kanoniczna**
macierz `Psi x Psi` albo ze zamknelo `QW-2191`.

`C15` robi cos wezszejszego:
- sprawdza, czy po przyjeciu control transport schema z `C14`
  da sie juz zapisac jawny formalny pullback packet
  `M_control = T_control^T H_PsiPsi T_control`,
- i czy repo ma juz coefficient-filled artefakty dla `H_PsiPsi` i `M_control`
  w zadeklarowanym (declared) scope, bez claimu physical canonicalization.

## Polityka zrodel

### Strict-admissible support

1. `QW-2190`
   - real Fourier control vectors `c1,s1,c2,s2` in `R^12`.
2. `QW-2164`
   - canonical Hessian shape `[13,13]`,
   - representative `Psi`-sector seed `eta0`.
3. `QW-2166`
   - exhaustive canonical Hessian,
   - cross-check `eta6`,
   - `linear_operator_matrix_matches_canonical_hessian=True`.
4. `QW-2180`
   - exact operator/Hessian identification on canonical carrier.
5. `C12`
   - minimal `Psi`-sector block extraction packet.
6. `C14`
   - control transport schema `mode basis -> Psi basis`.
7. `A10`
   - anti-overclaim boundary.

## Pytanie audytowe

Czy wolno juz powiedziec, ze:
- nawet bez physical canonicalization strict core daje juz formalny target
  assembled submatrix w trybie control-only,
- bo control transport schema z `C14` oraz canonical `Psi` Hessian carrier
  pozwalaja zapisac pullback
  `M_control = T_control^T H_PsiPsi T_control`,
- a wiec blocker nie dotyczy juz samej formuly assembly,
  tylko coefficient-filled canonical `H_PsiPsi`
  oraz restriction do orientation slice?

## Co daje `C14`

`C14` ustala:
- control transport schema `T_control : mode basis -> Psi basis` jest jawny
  na poziomie control identification,
- z control index-sets:
  - `I_mode_1 = {c1,s1}`
  - `I_mode_2 = {c2,s2}`.

To oznacza, ze mozna wybrac control basis:

- `B_control = (c1, s1, c2, s2)`

i traktowac jej kolumny jako wspolczynniki w canonical carrierze
`(psi0,...,psi11)`.

## Co daja `QW-2164/2166/2180`

Te bramki daja:
- canonical Hessian carrier dla `12` pol `Psi`,
- exhaustive operator/Hessian schema,
- exact identification operatora z Hessianem.

Najslabszy uczciwy obiekt kanoniczny brzmi:

- `H_PsiPsi`

jako formalny canonical `12 x 12` `Psi x Psi` block
wewnatrz exhaustive Hessian carrier.

Update (2026-03-16): repo eksportuje teraz jawny coefficient-filled canonical
blok `H_PsiPsi` na pelnym zadeklarowanym support-cie transportu (`R12`) oraz
wynikowy coefficient-filled declared control pullback `M_control` (`P476`).
To nadal jest ponizej physical canonicalization i ponizej host matching.

## Formalny packet po `C15`

Po zlozeniu `C12 + C14` mozna zapisac formalnie:

- `T_control in R^(12x4)`
- `H_PsiPsi in R^(12x12)`
- `M_control := T_control^T H_PsiPsi T_control in R^(4x4)`

To jest juz assembled target w trybie control-only.

Mozna go rozbic blokowo jako:

- `M_control = [[M_11, M_12], [M_21, M_22]]`

gdzie:
- `M_11` odpowiada parze `(c1,s1)`,
- `M_22` odpowiada parze `(c2,s2)`,
- `M_12`, `M_21` odpowiadaja sprzezeniom miedzy parami.

## Wynik `C15`

`C15` ustala:
- blocker `C14_B2` nie dotyczy juz samej formuly assembly w trybie control-only,
- strict core daje juz formalny pullback packet:
  - canonical `Psi x Psi` host block `H_PsiPsi`,
  - control transport `T_control`,
  - assembled control-only target `M_control = T_control^T H_PsiPsi T_control`.

## Redukcja frontu po `C15` (stan: 2026-03-16)

Po `C14` mielismy:

- `C14_B2 := no_assembled_Psi_x_Psi_submatrix_after_adopting_the_control_transport_schema`

Po `C15` (po update) najuczciwiej zapisac brak weziej jako:

- `C15_B2 := no_explicit_restriction_from_M_control_to_the_candidate_orientation_slice`

To jest realny postep redukcyjny:
- control-only assembly formula jest juz jawna,
- otwarty brak dotyczy coefficient filling i orientation restriction.

## Macierz wyniku

| Pytanie | Status po C15 | Uwagi |
|---|---|---|
| control-only assembly formula exists | `present_partial` | `M_control = T_control^T H_PsiPsi T_control` |
| coefficient-filled canonical `H_PsiPsi` exists | `yes` | `R12` (full declared support) |
| coefficient-filled `M_control` exists | `yes` | `P476` (declared control pullback) |
| restriction `M_control -> orientation slice` exists | `not_shown` | nadal brak |
| discharge of `C14_B2` | `reduced_not_closed` | blocker tylko zawężony |

## Czego `C15` nie ustala

`C15` nie ustala:
- ze export `M_control` jest juz selector-relevant (physical canonicalization transportu nadal otwarta),
- ze istnieje juz orientation-slice restriction,
- ze `C14_B2` ma PASS,
- ze `C9_B2` ma PASS.

## Anti-overclaim

`C15` nie twierdzi, ze:
- formalny pullback packet jest juz theorem-level block matchingiem,
- control-only assembly wystarcza do positivity descent,
- physical canonicalization transportu jest juz zalatwiona,
- uniqueness closure jest blisko.

## Produkt etapu

- pietnasty krok trzeciego mikrocyklu,
- jawny formalny packet assembly w trybie control-only,
- zawężenie `C14_B2` do braku coefficient-filled `H_PsiPsi`
  oraz braku restriction do orientation slice,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `C26`:
- zrobic pierwsze uczciwe zawężenie `M_control -> orientation slice`,
  bez claimu vanishing/cancellation i bez claimu host matching.
