# P1330 — Boundary-layer elimination attempt packet (PL)

Status: `BOUNDARY_LAYER_NOT_ELIMINATED`
As of: `2026-05-12`
Depends on: `P1329`

## Cel
Wykonać następny krok po `P1329`: zaatakować region graniczny
`|score_v4| < 0.03` i ocenić, czy można go wyeliminować miarowo jako
niedopuszczalny wkład residual loophole.

## Artefakt wykonawczy
- skrypt: `p1330_boundary_layer_elimination_attempt.py`
- raport: `generated/p1330_boundary_layer_elimination_attempt_report_v1.json`

## Wynik
- `boundary_ratio = 0.1983` (~19.8%),
- `exclusion_by_measure_supported = false`,
- status: `BOUNDARY_TOO_LARGE_FOR_EXCLUSION`.

## Decyzja profesorska
Próba eliminacji warstwy granicznej przez argument "małej miary" nie działa na
obecnym modelu: region graniczny jest zbyt duży, aby uczciwie go odrzucić.

## Konsekwencja
- Globalne L2 nadal `NOT_EXPORTED`.
- `QW-2191` strict nadal `NOT_CLOSED`.
- Potrzebny inny mechanizm: nie miarowy, lecz strukturalny (np. dodatkowe
  constraints/regularity, które rozbiją region graniczny).

## Czego dokument nie twierdzi
- Nie twierdzi globalnego L2 discharge.
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.

## Rekomendowany następny uczciwy krok
Uruchomić **P1331 structural boundary resolution attempt**:
1. dodać jawny warunek regularności strict dla near-boundary,
2. sprawdzić, czy warunek redukuje region graniczny bez sztucznego cięcia,
3. jeśli tak, ponowić próbę globalnego L2 export.

## Dla laika
Sprawdziliśmy, czy problematyczny "pas graniczny" jest na tyle mały, by go
uczciwie pominąć. Nie jest — jest za duży. Trzeba więc znaleźć głębszą regułę,
która porządkuje ten pas, zamiast go ignorować.
