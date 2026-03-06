# C19 Generator-Level All-Rows Source Audit

Status: `C19_EXECUTED_GENERATOR_LEVEL_ALL_ROWS_SOURCE_AUDIT_NO_FALSE_PASS`
As of: `2026-03-06`

## Cel

Po `C18` najwezszy aktywny blocker brzmial:

- `C18_B1 := no_explicit_serialized_12_row_export_table_for_the_Psi_family_despite_the_existing_finite_family_witness_packet`

`C19` nie udaje, ze strict core ma juz jawna, trwala tabele
`12` rows `Psi`.

`C19` robi cos wezszejszego:
- sprawdza, czy strict core ma juz **generator-level all-rows source**
  dla calej rodziny `Psi`,
- i czy blocker da sie zawezic do braku trwalego, jawnie
  zmaterializowanego artefaktu serializacji,
  a nie do braku samego zrodla dla wszystkich rows.

## Polityka zrodel

### Strict-admissible support

1. `QW-2165`
   - generator liczy `eom_psi[i]` dla wszystkich `12` fields `Psi`,
   - report przechowuje:
     - `lagrangian_density`,
     - trzy sample rows:
       - `sample_eom_psi0`,
       - `sample_eom_psi6`,
       - `sample_eom_psi11`,
   - oraz exhaustive flags dla calej rodziny.
2. `QW-2166`
   - generator buduje pelny Hessian dla wszystkich `13` pol,
   - report przechowuje:
     - `potential_total`,
     - `hessian_shape = [13, 13]`,
     - sample linearized equations,
   - oraz exhaustive Hessian/operator flags.
3. `QW-2180`
   - exact operator/Hessian identification.
4. `C18`
   - finite family-level witness packet.
5. `A10`
   - anti-overclaim boundary.

## Pytanie audytowe

Czy wolno juz powiedziec, ze:
- strict core ma generator-level source dla wszystkich `12` rows `Psi`,
- ale nadal nie ma jawnie zmaterializowanego `12`-row export artifact,
- a wiec blocker dotyczy juz tylko warstwy persisted serialization,
  nie braku samego all-rows source?

## Co daje `QW-2165`

`QW-2165` daje dwie rzeczy naraz:

1. w pliku generatora:
   - `eom_psi = [euler_lagrange(l_density, psi[i]) for i in range(N)]`
   - czyli generator-level source dla wszystkich rows `Psi`;
2. w raporcie:
   - `n_psi_fields = 12`,
   - pelne `lagrangian_density`,
   - tylko trzy sample rows:
     - `psi0`,
     - `psi6`,
     - `psi11`.

To jest dokladnie roznica, ktora `C19` chce nazwac:
- source juz jest generator-level exhaustive,
- persisted export nadal jest tylko probkowany.

## Co daje `QW-2166`

`QW-2166` daje warstwe analogiczna po stronie Hessian/operator:

- w pliku generatora:
  - `hessian = sp.hessian(potential_total, fields)`,
- w raporcie:
  - `n_psi_fields = 12`,
  - `potential_total`,
  - `hessian_shape = [13, 13]`,
  - exhaustive flags:
    - `hessian_constructed_for_all_13_fields=True`,
    - `linearized_eom_executed_for_all_13_fluctuation_fields=True`.

To wzmacnia wniosek:
- brak nie dotyczy juz samej obecnosci generator-level source,
- tylko jego jawnej, finalnej serializacji w artefakcie eksportowym.

## Wynik `C19`

`C19` ustala:
- strict core ma juz generator-level all-rows source
  dla calej rodziny `12` rows `Psi`,
- finite witness packet z `C18` nie byl tylko slaba heurystyka,
  bo pod spodem istnieje juz exhaustive source computation,
- aktualny blocker nie dotyczy juz braku source,
  tylko braku persisted serialization artifact.

## Redukcja frontu po `C19`

Po `C18` mielismy:

- `C18_B1 := no_explicit_serialized_12_row_export_table_for_the_Psi_family_despite_the_existing_finite_family_witness_packet`
- `C18_B2 := no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice`

Po `C19` najuczciwiej zapisac to weziej jako:

- `C19_B1 := no_explicit_persisted_12_row_serialization_artifact_even_though_generator_level_all_rows_source_is_present`
- `C19_B2 := no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice`

To jest realny postep redukcyjny:
- generator-level exhaustive source jest juz obecny,
- nadal brak finalnego, jawnego artefaktu eksportowego.

## Macierz wyniku

| Pytanie | Status po C19 | Uwagi |
|---|---|---|
| generator-level all-rows source exists | `present_partial` | `QW-2165` i `QW-2166` licza obiekty exhaustive |
| persisted serialized 12-row export exists | `not_shown` | nadal brak |
| persisted full `12 x 12` table exists | `not_shown` | nadal brak |
| restriction to candidate orientation slice exists | `not_shown` | nadal brak |
| discharge of `C18_B1` | `reduced_not_closed` | blocker tylko zawężony |

## Czego `C19` nie ustala

`C19` nie ustala:
- ze wszystkie `12` rows sa juz jawnie wypisane w jednym artefakcie,
- ze istnieje juz finalna tabela `12 x 12`,
- ze control pullback zostal juz policzony liczbowo lub symbolicznie do konca,
- ze restriction do candidate orientation slice istnieje,
- ze `C18_B1` ma PASS,
- ze `C18_B2` ma PASS.

## Anti-overclaim

`C19` nie twierdzi, ze:
- generator-level source jest rownowazny persisted export artifact,
- obecny report rowna sie gotowej row-by-row serializacji,
- z samej obecnosci `lagrangian_density` i `potential_total`
  wynika juz jawna tabela `12` rows.

## Produkt etapu

- dziewietnasty krok trzeciego mikrocyklu,
- generator-level all-rows source audit dla rodziny `Psi`,
- zawężenie `C18_B1` do braku persisted serialization artifact,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `C20`:
- sprawdzic, czy strict core pozwala juz na jawny
  materialized `12`-row serialization packet z istniejacego
  generator-level source,
- albo jawnie potwierdzic, ze taki packet nadal nie jest obecny.
