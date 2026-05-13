# P1475 — S4.25 QW-2191 SP1 Window Stress Replay (PL)

Status: `P1475_EXECUTED_QW2191_SP1_WINDOW_STRESS_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Gęsty replay przy dolnej granicy okna SP1 z P1474, aby sprawdzić numeryczną
stabilność granicy i uniknąć przypadkowego overclaim.

## Rygor

- bez legacy bridge,
- bez strict-core closure claim,
- tylko local-only boundary stress evidence.
