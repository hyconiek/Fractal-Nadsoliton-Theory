# P1329 — Residual-loophole elimination proof attempt (PL)

Status: `PARTIAL_CONDITIONAL_SUPPORT_NOT_GLOBAL_EXPORT`
As of: `2026-05-12`
Depends on: `P1328`

## Cel
Wykonać rekomendowany krok `P1328`: formalny case split dla `z2/eps` i ocenić,
czy da się wyeksportować lemę L2 (`Residual_loophole_elimination`).

## Artefakt wykonawczy
- skrypt: `p1329_residual_loophole_elimination_proof_attempt.py`
- raport: `generated/p1329_residual_loophole_elimination_proof_attempt_report_v1.json`

## Metoda
- siatka `phase x amp`,
- residual shift `0.015 * eps * z2`,
- domena admissible: `|score_v4| >= margin_delta` z `margin_delta = 0.03`,
- test: czy residual może odwrócić znak kierunku.

## Wynik
- `admissible_points = 650`,
- `violations = 0` w domenie admissible,
- `l2_proven_on_admissible_margin_domain = true`,
- `global_l2_exportable = false`.

## Decyzja profesorska
To jest **warunkowe wsparcie** L2, ale nie pełny eksport globalny.

Mamy dowód braku odwrócenia znaku na domenie z marginesem `|score_v4|>=0.03`,
ale nie mamy jeszcze globalnego twierdzenia dla całej przestrzeni (szczególnie
near-boundary region).

## Konsekwencja
- `QW-2191` strict nadal `NOT_CLOSED`.
- L2 można oznaczyć jako `CONDITIONAL_SUPPORT`, nie `FORMALLY_EXPORTED`.

## Czego dokument nie twierdzi
- Nie twierdzi globalnego L2 discharge.
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.

## Rekomendowany następny uczciwy krok
Uruchomić **P1330 boundary-layer elimination attempt**:
1. zaatakować near-boundary region `|score_v4| < 0.03`,
2. albo wykazać jego niedopuszczalność strict-core,
3. albo domknąć globalne ograniczenie residual shift.

## Dla laika
Udało się pokazać, że w "bezpiecznej strefie" kompasu pokrętło `z2/eps` nie
potrafi zmienić kierunku. Ale w cienkiej strefie granicznej wciąż nie mamy
pełnej gwarancji, więc nie można jeszcze ogłosić pełnego domknięcia.
