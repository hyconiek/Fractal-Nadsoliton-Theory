# P1170 Neighbor Candidate E2E Screen Packet

Status: `P1170_EXECUTED_NEIGHBOR_CANDIDATE_E2E_SCREEN_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Po `P1169` wykonujemy krok porównawczy: sprawdzić stabilność wyboru top
kandydata na małym sąsiedztwie `(sigma,kappa)` i potwierdzić przejście E2E.

## Professor-level decision

Buduję 5 kandydatów sąsiednich wokół punktu top z `P1166` i dla każdego:

1. liczę metryki proxy (`positive/negative/sign_change`),
2. uruchamiam pełny `P1151` pipeline,
3. tworzę ranking mieszany (pipeline pass + burden metryk).

## Artifacts

- script:
  `p1170_neighbor_candidate_e2e_screen.py`
- summary:
  `generated/p1170_neighbor_candidate_e2e_screen_summary.json`
- generated neighbor candidates:
  `generated/p1170_neighbor_*.json`

## Result

Przesiano 5 kandydatów; summary zawiera pełne porównanie i `top_recommendation`.

## Honest boundary

`P1170` to walidacja i porównanie kandydatów roboczych, bez claimu closure i
bez `QW-2191` discharge.
