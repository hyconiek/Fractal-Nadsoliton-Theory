# P1167 Stability Extension Probe Packet

Status: `P1167_EXECUTED_STABILITY_EXTENSION_PROBE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać następny uczciwy krok po `P1166`: sprawdzić stabilność lokalnego
`CANDIDATE_PASS_ONLY` poza domeną bazową i przy małych perturbacjach parametrów.

## Professor-level decision

Dla najlepszego punktu z `P1166`:

1. test bazowy na `[0,24]`,
2. test rozszerzony na `[0,48]`,
3. test odporności na perturbacje `(domega,dphi,dsigma) in {-0.01,0,0.01}^3`
   (27 przypadków) na `[0,48]`.

## Result

- rozszerzony test `[0,48]` zachowuje `sign_change_count=0` i `negative_count=0`,
- odporność perturbacyjna: `robust_cases = 27/27` (`robust_fraction=1.0`).

To wzmacnia status lokalnego kandydata jako stabilnego **proxy-pass**, nadal bez
claimu closure.

## Artifacts

- script:
  `p1167_stability_extension_probe.py`
- summary:
  `generated/p1167_stability_extension_probe_summary.json`

## Honest boundary

`P1167` to test stabilności proxy, nie dowód domknięcia teorii ani discharge
`QW-2191`.
