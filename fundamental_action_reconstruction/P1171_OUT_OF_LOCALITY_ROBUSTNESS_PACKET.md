# P1171 Out-of-Locality Robustness Packet

Status: `P1171_EXECUTED_OUT_OF_LOCALITY_ROBUSTNESS_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Po `P1170` wykonujemy test poza-lokalny: czy top kandydat pozostaje stabilny
po rozszerzeniu domeny i większych perturbacjach parametrów.

## Professor-level decision

Dla top kandydata `P1170`:

1. test bazowy na `[0,72]`,
2. perturbacje 4D:
   - `domega in {-0.02,0,0.02}`
   - `dphi   in {-0.03,0,0.03}`
   - `dsigma in {-0.2,0,0.2}`
   - `dkappa in {-0.04,0,0.04}`

(81 przypadków).

## Result

- baza `[0,72]`: `sign_change_count=0`, `negative_count=0`,
- odporność poza-lokalna: `robust_cases = 48/81` (`robust_fraction ≈ 0.593`).

Wniosek: kandydat jest stabilny w rdzeniu, ale nie globalnie odporny na cały
zakres perturbacji — potrzebna dalsza kalibracja regionu roboczego.

## Artifacts

- script:
  `p1171_out_of_locality_robustness_probe.py`
- summary:
  `generated/p1171_out_of_locality_robustness_probe_summary.json`

## Honest boundary

`P1171` to test odporności proxy, nie dowód closure ani `QW-2191` discharge.
