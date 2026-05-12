# P1336 — Internal-source v2 tie-break family search packet (PL)

Status: `V2_EXPORTABLE_WITH_DEGENERATE_ALPHA0_WARNING`
As of: `2026-05-12`
Depends on: `P1335`

## Cel
Wykonać rekomendowany krok z `P1335`: przeszukać rodzinę tie-break i znaleźć
wariant strict-side zgodny z operacyjną regułą branch-admissibility.

## Artefakt wykonawczy
- skrypt: `p1336_internal_source_v2_tiebreak_family_search.py`
- raport: `generated/p1336_internal_source_v2_tiebreak_family_search_report_v1.json`

## Wynik
- najlepszy parametr: `alpha = 0.0`,
- `mismatches = 0`,
- `strict_internal_source_exportable = true`,
- status: `V2_EXPORTABLE`.

## Decyzja profesorska
Mamy pierwszy strict-side kompatybilny kandydat źródła gałęzi.

Uwaga metodologiczna: optimum `alpha=0` oznacza neutralny tie-break (bez
dodatkowego przesunięcia near-boundary). To jest zgodne, ale wymaga dalszej
interpretacji theorem-level, aby uniknąć zarzutu "trywialnej zgodności".

## Konsekwencja
- Bariera `strict_internal_source_for_A_branch_v1_exported` przechodzi na `true`.
- Nadal brak globalnego formal export L2, więc `QW-2191 strict` pozostaje
  `NOT_CLOSED` do czasu pełnego O3 checker pass.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że sam search parametryczny kończy sprawę theorem-level.

## Rekomendowany następny uczciwy krok
Uruchomić **P1337 full O3 checker refresh**:
1. zaktualizować obligation checker o nowy export internal-source,
2. sprawdzić, czy jedyną brakującą pozycją pozostaje globalne L2,
3. jeśli tak, skoncentrować się na finalnym L2 global export packet.

## Dla laika
W drugiej próbie znaleźliśmy ustawienie, które nie kłóci się z działającą regułą.
To ważny postęp. Ale nadal brakuje ostatniego formalnego "stempla" globalnego,
więc pełnego domknięcia strict jeszcze nie ogłaszamy.
