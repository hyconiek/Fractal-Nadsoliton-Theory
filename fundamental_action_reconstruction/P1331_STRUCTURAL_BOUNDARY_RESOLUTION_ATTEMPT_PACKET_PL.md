# P1331 — Structural boundary resolution attempt packet (PL)

Status: `STRUCTURAL_RESOLUTION_WEAK_NOT_SUFFICIENT`
As of: `2026-05-12`
Depends on: `P1330`

## Cel
Po nieudanej eliminacji miarowej (`P1330`) sprawdzić, czy warunek regularności
(gradient floor) redukuje problematyczny region graniczny.

## Artefakt wykonawczy
- skrypt: `p1331_structural_boundary_resolution_attempt.py`
- raport: `generated/p1331_structural_boundary_resolution_attempt_report_v1.json`

## Metoda
- warunek admissibility: `||grad score_v4|| >= gmin`, `gmin = 0.35`,
- granica: `|score_v4| < delta`, `delta = 0.03`,
- ocena: udział granicy w domenie admissible.

## Wynik
- `boundary_ratio_in_admissible = 0.1983` (~19.8%),
- `structural_resolution_supported = false`,
- status: `STRUCTURAL_SUPPORT_WEAK`.

## Decyzja profesorska
Aktualny warunek regularności nie rozwiązuje problemu granicznego — praktycznie
nie redukuje udziału boundary-layer do poziomu akceptowalnego dla globalnego L2.

## Konsekwencja
- Globalne `Residual_loophole_elimination` nadal `NOT_EXPORTED`.
- `QW-2191` strict nadal `NOT_CLOSED`.
- Potrzebny bardziej selektywny warunek strukturalny (np. warunek krzywizny,
  nie tylko gradient floor).

## Czego dokument nie twierdzi
- Nie twierdzi globalnego L2 discharge.
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.

## Rekomendowany następny uczciwy krok
Uruchomić **P1332 curvature-aware boundary suppression attempt**:
1. dodać warunek krzywizny/second-order stability,
2. sprawdzić redukcję boundary-layer,
3. jeśli skuteczne, zaktualizować L2 export readiness.

## Dla laika
Dodaliśmy nową regułę "gładkości" kompasu, ale ona nie pomogła — problematyczny
pas graniczny nadal jest duży. Trzeba użyć bardziej precyzyjnej reguły, która
patrzy nie tylko na nachylenie, ale też na krzywiznę.
