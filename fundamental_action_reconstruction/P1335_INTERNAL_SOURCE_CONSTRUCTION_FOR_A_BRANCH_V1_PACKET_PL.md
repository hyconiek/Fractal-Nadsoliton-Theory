# P1335 — Internal-source construction for `A_branch_v1` packet (PL)

Status: `INTERNAL_SOURCE_ATTEMPT_INCOMPATIBLE`
As of: `2026-05-12`
Depends on: `P1334`

## Cel
Wykonać rekomendowany krok z `P1334`: skonstruować strict-side kandydat
mechanizmu wyboru gałęzi near-boundary i sprawdzić jego kompatybilność z
operacyjnym zachowaniem `A_branch_v1`.

## Artefakt wykonawczy
- skrypt: `p1335_internal_source_construction_for_a_branch_v1.py`
- raport: `generated/p1335_internal_source_construction_for_a_branch_v1_report_v1.json`

## Wynik
- `boundary_points = 617`,
- `mismatches_vs_axiom = 205`,
- `candidate_internal_source_compatible_with_axiom = false`,
- status: `INTERNAL_SOURCE_CANDIDATE_INCOMPATIBLE`.

## Decyzja profesorska
Pierwsza konstrukcja strict-side źródła nie przechodzi: zbyt często nie zgadza
się z regułą branch-admissibility używaną operacyjnie w near-boundary.

## Konsekwencja
- `strict_internal_source_exportable = false` na obecnej wersji.
- `QW-2191 strict` pozostaje `NOT_CLOSED`.
- Potrzebna nowa konstrukcja (`v2`) mechanizmu strict-side.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że strict derivation jest niemożliwa.

## Rekomendowany następny uczciwy krok
Uruchomić **P1336 internal-source v2 with constrained tie-break family**:
1. zdefiniować rodzinę tie-break opartą na lokalnych inwariantach,
2. przeszukać parametry pod minimalizację mismatch,
3. dopuścić eksport dopiero przy `mismatches = 0` i replay pass.

## Dla laika
Próbowaliśmy zbudować "wewnętrzny automat" teorii, który sam wybiera gałąź.
Ta pierwsza wersja często nie zgadza się z dotychczas działającą regułą,
więc nie nadaje się jeszcze do ogłoszenia pełnego domknięcia strict.
