# P1332 — Curvature-aware boundary suppression attempt packet (PL)

Status: `CURVATURE_SUPPORT_WEAK_NOT_SUFFICIENT`
As of: `2026-05-12`
Depends on: `P1331`

## Cel
Wykonać krok `P1332`: sprawdzić, czy warunek drugiego rzędu (krzywizna)
redukuje boundary-layer skuteczniej niż gradient-only z `P1331`.

## Artefakt wykonawczy
- skrypt: `p1332_curvature_aware_boundary_suppression_attempt.py`
- raport: `generated/p1332_curvature_aware_boundary_suppression_attempt_report_v1.json`

## Metoda
- score korygowany o termin krzywizny:
  `|score_v4| - 0.5*|d2 score/d phase^2|*h^2`, gdzie `h=0.02`,
- boundary w domenie admissible: skorygowany score `< delta`, `delta=0.03`.

## Wynik
- `boundary_ratio_in_admissible = 0.1988` (~19.9%),
- `suppression_supported = false`,
- status: `CURVATURE_SUPPORT_WEAK`.

## Decyzja profesorska
Korekta krzywiznowa w tej postaci nie rozwiązuje problemu boundary-layer.
Wynik praktycznie nie poprawia stanu względem `P1331`.

## Konsekwencja
- Globalne L2 nadal `NOT_EXPORTED`.
- `QW-2191` strict nadal `NOT_CLOSED`.
- Trzeba przejść z lokalnych korekt geometrycznych do mocniejszego mechanizmu
  logiczno-strukturalnego (np. jawna reguła dopuszczalności gałęzi).

## Czego dokument nie twierdzi
- Nie twierdzi globalnego L2 discharge.
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.

## Rekomendowany następny uczciwy krok
Uruchomić **P1333 branch-admissibility axiom test (non-strict tagged)**:
1. jawnie wprowadzić kandydat reguły dopuszczalności gałęzi near-boundary,
2. przetestować wpływ na globalne L2,
3. utrzymać etykietę non-strict, dopóki reguła nie będzie wyprowadzona strict-side.

## Dla laika
Spróbowaliśmy bardziej zaawansowanej poprawki (z krzywizną), ale problematyczny
pas graniczny prawie się nie zmienił. To znaczy, że same sztuczki geometryczne
nie wystarczą — trzeba dodać mocniejszą zasadę wyboru gałęzi.
