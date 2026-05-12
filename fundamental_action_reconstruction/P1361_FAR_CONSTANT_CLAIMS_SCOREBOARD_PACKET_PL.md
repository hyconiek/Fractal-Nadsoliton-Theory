# P1361 FAR Constant Claims Scoreboard Packet (PL)

Status: `P1361_EXECUTED_SCOREBOARD_NO_FALSE_PASS`
As of: `2026-05-12`
Artifacts:
- `generated/p1361_far_constant_claims_scoreboard_summary.json`

## Cel

Wykonać profesorski następny krok po `P1360`: zamiast ogólnej narracji zrobić tablicę wyników (scoreboard) wszystkich głównych claimów stałych fizycznych w FAR.

## Co zrobiono

Uruchomiono `p1361_far_constant_claims_scoreboard_checkpoint.py`, który:

1. klasyfikuje claimy do: `strict_verified / strict_candidate / nonstrict_proxy / legacy_only`,
2. sprawdza obecność artefaktów dowodowych,
3. zapisuje licznik statusów i next priority.

## Wynik bieżący

Aktualny licznik:

- `strict_verified = 0`
- `strict_candidate = 3`
- `nonstrict_proxy = 1`
- `legacy_only = 1`

To jest uczciwy obraz stanu: mamy silne kandydaty, ale brak jeszcze pełnego strict-verified closure dla stałych fizycznych.

## Decyzja profesorska

Następny uczciwy krok: `P1362`

1. dla każdego `strict_candidate` dołożyć obowiązkowy benchmark residuali,
2. dołożyć jawny budget niepewności,
3. podnieść status na `strict_verified` tylko przy residual pass + stability pass.

## Dla laika

Zrobiliśmy „tablicę wyników” teorii: co już jest naprawdę potwierdzone, a co jeszcze kandydatem.
Na dziś kandydatów jest dużo, ale pełnych potwierdzeń jeszcze nie ma.
To jest dobra nauka, bo pokazuje prawdę, a nie marketing.
