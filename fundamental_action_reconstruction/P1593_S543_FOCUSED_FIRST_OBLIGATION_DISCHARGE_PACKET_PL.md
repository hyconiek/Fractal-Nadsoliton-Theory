# P1593 / S543 Focused First Obligation Discharge Packet (PL)

Status: `P1593_EXECUTED_FOCUSED_DISCHARGE_STEP`
As of: `2026-05-14`

## Cel

Wykonać kolejny uczciwy krok po syntezie progressu:

1. pobrać pierwszy brak z listy closure,
2. zaatakować go formalnym kandydatem discharge,
3. utrzymać pełny strict tor i status `OPEN` gdy warunki nie są spełnione.

## Wynik

- Eksportuje focused-obligation status i metryki dowodowe.
- Przesuwa listę braków tylko gdy warunki candidate-pass są spełnione.

## Artefakt

- `generated/p1593_s543_focused_first_obligation_discharge_summary.json`

## Następny uczciwy krok

`P1594`: discharge kolejnego brakującego obiektu (typowo global stability theorem object).
