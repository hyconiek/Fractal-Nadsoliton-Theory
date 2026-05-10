# P1166 Sigma-Kappa Sweep Packet

Status: `P1166_EXECUTED_SIGMA_KAPPA_SWEEP_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać zalecany krok po `P1165`: 2D skan `(sigma, kappa)` dla klasy
asymetrycznego członu selekcyjnego.

## Professor-level decision

Przy stałych parametrach bazowych `(omega,phi,beta,eta)` skanuję siatkę:

- `sigma in {0.2,0.4,0.6,0.8,1.0,1.2}`
- `kappa in {0.04,0.08,0.12,0.16,0.20,0.30}`

(36 punktów).

## Result

- `grid_size = 36`,
- `zero_sign_change_count = 4`.

To jest pierwszy lokalny sygnał, że część konfiguracji może osiągać
`sign_change_count=0` na domenie `[0,24]`.

Jednocześnie brak claimu closure: trzeba jeszcze sprawdzić pełne kryteria
fizyczne i stabilność poza lokalnym skanem.

## Artifacts

- script:
  `p1166_sigma_kappa_sweep.py`
- summary:
  `generated/p1166_sigma_kappa_sweep_summary.json`

## Honest boundary

`P1166` to wynik lokalnego skanu kandydatów, nie dowód domknięcia teorii i nie
`QW-2191` discharge.
