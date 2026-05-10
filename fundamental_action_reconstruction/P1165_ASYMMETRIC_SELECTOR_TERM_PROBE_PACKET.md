# P1165 Asymmetric Selector Term Probe Packet

Status: `P1165_EXECUTED_ASYMMETRIC_SELECTOR_TERM_PROBE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Po lokalnych negatywnych wynikach faza+tłumienie (`P1163`,`P1164`) wykonujemy
następny uczciwy krok: test nowej klasy z jawnym asymetrycznym składnikiem
selekcyjnym strict-side.

## Professor-level decision

Badana klasa:

```text
N_sigma(d) = cos(omega*d+phi) + sigma*(1-exp(-kappa*d)), sigma>=0
K_sigma(d) = exp(-alpha*d) * N_sigma(d) / (1+beta*d^eta)
```

Przy stałych parametrach bazowych z `P1163-best` skanujemy `sigma`.

## Result

W badanym zakresie `sigma`:

- `zero_sign_change_count = 0`,
- najlepszy przypadek nadal proxy `BLOCKED`.

To oznacza, że w tej lokalnej wersji asymetrycznego członu nie uzyskano
przejścia do `sign_change_count=0`.

## Artifacts

- script:
  `p1165_asymmetric_selector_term_probe.py`
- summary:
  `generated/p1165_asymmetric_selector_term_probe_summary.json`

## Honest boundary

`P1165` to test klasy kandydatów, nie dowód domknięcia teorii.
