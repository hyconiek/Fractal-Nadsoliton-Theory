# P1161 Physical Candidate Sweep Packet

Status: `P1161_EXECUTED_PHYSICAL_CANDIDATE_SWEEP_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Po kontrolnych testach `P1160` wykonujemy bardziej fizycznie sensowny krok:
porównanie kandydatów strict-side z jawnie zadanymi hintami parametrów
`omega/phi/beta/eta` bez sztucznego `force_fail_stage`.

## Professor-level decision

Dodaję lokalny sweep kandydatów (`P1161`):

1. dwa kandydaty fizyczne A/B,
2. metryki typu `P1146` (negative_count, sign_change_count),
3. ranking po najmniejszym burden oscylacyjnym.

## Artifacts

- candidates:
  - `generated/p1161_candidate_physical_a.json`
  - `generated/p1161_candidate_physical_b.json`
- registry:
  - `generated/p1161_candidate_registry_physical.json`
- sweep script:
  - `p1161_candidate_probe_sweep.py`
- summary:
  - `generated/p1161_candidate_probe_sweep_summary.json`

## Result

Skrypt wylicza porównanie kandydatów i wskazuje `top_recommendation` na
podstawie jawnych kryteriów metrycznych.

## Honest boundary

`P1161` to nadal etap badania kandydatów (metodologia), bez claimu closure i
bez discharge `QW-2191`.
