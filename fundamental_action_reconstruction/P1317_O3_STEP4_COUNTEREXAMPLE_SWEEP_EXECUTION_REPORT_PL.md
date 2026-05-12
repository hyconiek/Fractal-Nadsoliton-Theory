# P1317 — O3.4 Counterexample sweep execution report (strict-only, PL)

Status: `O3_4_EXECUTED_COUNTEREXAMPLE_SWEEP_PARTIAL`
As of: `2026-05-12`
Depends on: `P1316`

## Cel
Wykonać O3.4 operacyjnie (nie tylko plan): uruchomić reproducowalny sweep
kontrprzykładów na parach C1–C4 i sprawdzić, czy pojawia się legalna kontr-gałąź
kierunku w strict lane.

## Artefakt wykonawczy
- skrypt: `p1317_o34_counterexample_sweep.py`
- raport JSON: `generated/p1317_o34_counterexample_sweep_report_v1.json`

## Wynik uruchomienia
- `o3_4_status = PASS_NO_COUNTEREXAMPLE`
- `COUNTEREXAMPLE_FOUND = 0`
- `RESIDUAL_AMBIGUITY > 0`
- `nonuniqueness_residual_status = OPEN`
- `qw2191_strict_status = NOT_CLOSED`

Interpretacja profesorska:
1. nie znaleziono twardego kontrprzykładu łamiącego bieżące tolerancje,
2. ale nadal istnieje nierozstrzygnięta wieloznaczność residualna,
3. więc uczciwie nie wolno ogłaszać domknięcia `QW-2191`.

## Konsekwencja logiczna dla O3
- O3.4 usuwa część ryzyka (brak jawnego kontrprzykładu),
- ale nie domyka O3, bo O3 wymaga także wyzerowania residual ambiguity.

## Czego dokument nie twierdzi
- Nie twierdzi `QW-2191 = CLOSED`.
- Nie twierdzi strict-core selector closure.
- Nie twierdzi globalnego ToE closure.

## Rekomendowany następny uczciwy krok
Wykonać **O3.5**: niezależny replay + adversarial audit (inna kolejność,
zmienne seed/control), ze szczególnym naciskiem na gałąź `open(Z2/eps)` z C3.

## Dla laika
To tak, jakbyśmy testowali różne kompasy i nie znaleźli żadnego oczywiście
zepsutego. Ale część kompasów dalej ma "luz" w mechanizmie, więc nie możemy
jeszcze powiedzieć, że wszystkie zawsze pokażą dokładnie to samo.
