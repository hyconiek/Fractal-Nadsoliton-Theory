# C29 Local Projector Formula Packet

Status: `C29_EXECUTED_LOCAL_PROJECTOR_FORMULA_PACKET_NO_FALSE_PASS`
As of: `2026-03-06`

## Cel

Po `C28` najwezszy aktywny blocker brzmi:

- `C28_B1 := no_explicit_serialized_local_orbit_frame_projector_formula_or_global_gluing_rule_for_the_control_coordinate_quotient_candidate`

`C29` nie probuje udawac, ze strict core ma juz globalny quotient map.

`C29` robi cos wezszejszego:
- sprawdza, czy w lokalnej parze `(c_ref,s_ref)` da sie juz jawnie zapisac
  formule projektora:
  - na kierunek orbit tangent,
  - oraz na kierunek poprzeczny po lokalnym quociencie.

## Polityka zrodel

### Strict-admissible support

1. `C4`
   - lokalna geometria orbity `O(2)`,
   - `u(theta) = cos(theta) c_ref + sin(theta) s_ref`.
2. `C28`
   - local orbit-frame schema:
     - tangent direction,
     - transverse mismatch direction.
3. `C14`
   - control basis osadzona w canonical carrierze.
4. `C15`
   - control-only pullback carrier `M_control`.
5. `A10`
   - anti-overclaim boundary.

## Pytanie audytowe

Czy wolno juz powiedziec, ze strict core ma packet-ready
**jawna lokalna formule projektora** na parze `(c_ref,s_ref)`:

- `P_tan(theta)` na kierunek orbit tangent,
- `P_red(theta)` na lokalny reduced direction po quociencie,

nawet jesli nadal brak:
- globalnego sklejenia miedzy parami,
- finalnego wydobycia dwuwymiarowej orientation slice?

## Lokalna rama orbity

W lokalnej bazie `((c_ref),(s_ref))` zapisujemy:

- `e(theta) := (cos(theta), sin(theta))^T`,
- `tau(theta) := (-sin(theta), cos(theta))^T`.

Mamy:
- `e(theta)` i `tau(theta)` sa ortonormalne,
- `tau(theta)` jest styczny do orbity,
- `e(theta)` jest lokalnym kandydatem na kierunek poprzeczny
  po lokalnym quociencie.

## Jawne formuly projektorow

Lokalny projector na tangent direction:

```text
P_tan(theta) = tau(theta) tau(theta)^T
             = [[ sin^2(theta),   -sin(theta)cos(theta)],
                [ -sin(theta)cos(theta), cos^2(theta) ]]
```

Lokalny projector na reduced transverse direction:

```text
P_red(theta) = e(theta) e(theta)^T
             = [[ cos^2(theta),  cos(theta)sin(theta)],
                [ cos(theta)sin(theta), sin^2(theta) ]]
```

oraz:

```text
P_tan(theta) + P_red(theta) = I_2
P_tan(theta)^2 = P_tan(theta)
P_red(theta)^2 = P_red(theta)
P_tan(theta) P_red(theta) = 0
```

To jest juz jawna serialized local projector formula.

## Najmocniejszy uczciwy wniosek po `C29`

Po zlozeniu `C4 + C28 + C14 + C15` najuczciwiej zapisac:

- lokalny orbit-frame quotient candidate nie jest juz tylko schema,
- strict core ma juz packet-ready jawna formule lokalnych projektorow
  `P_tan(theta)` i `P_red(theta)`,
- ale nadal brak:
  - globalnego gluing rule miedzy parami `(c1,s1)` i `(c2,s2)`,
  - finalnego extraction map do dwuwymiarowej orientation slice.

## Redukcja frontu po `C29`

Po `C28` mielismy:

- `C28_B1 := no_explicit_serialized_local_orbit_frame_projector_formula_or_global_gluing_rule_for_the_control_coordinate_quotient_candidate`
- `C26_B2 := no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane`

Po `C29` najuczciwiej zapisac to weziej jako:

- `C29_B1 := no_explicit_pair_to_pair_global_gluing_rule_assembling_the_local_reduced_lines_into_a_single_reduced_control_plane`
- `C26_B2 := no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane`

To jest realny postep redukcyjny:
- lokalna serialized formula juz istnieje,
- otwarty pozostaje globalny gluing i finalna slice extraction.

## Macierz wyniku

| Pytanie | Status po C29 | Uwagi |
|---|---|---|
| explicit local tangent projector formula exists | `present_partial` | tak |
| explicit local reduced projector formula exists | `present_partial` | tak |
| global pair-to-pair gluing rule exists | `not_shown` | nadal brak |
| final basis-level slice extraction exists | `not_shown` | nadal brak |
| discharge of `C28_B1` | `reduced_not_closed` | tylko zawężenie |

## Czego `C29` nie ustala

`C29` nie ustala:
- ze lokalne projektory sa juz sklejone globalnie,
- ze istnieje globalny reduced control plane,
- ze istnieje finalna dwuwymiarowa orientation slice,
- ze `M_control` zostalo juz zredukowane do finalnego slice operatora,
- ze `C28_B1` ma PASS.

## Anti-overclaim

`C29` nie twierdzi, ze:
- lokalna formula projektora daje juz globalny quotient,
- selector track jest blisko theorem-level closure,
- finalny slice extraction zostal rozwiazany.

## Produkt etapu

- dwudziesty dziewiaty krok trzeciego mikrocyklu,
- jawna serialized local projector formula,
- zawężenie `C28_B1` do braku globalnego gluing rule,
- utrzymany brak falszywego PASS.

## Nastepny krok

Naturalnym kolejnym ruchem jest `C30`:
- sprawdzic, czy strict core ma juz packet-ready pair-to-pair gluing rule
  dla dwoch lokalnych reduced lines,
- albo jawnie potwierdzic, ze globalny reduced control plane nadal nie jest eksportowany.
