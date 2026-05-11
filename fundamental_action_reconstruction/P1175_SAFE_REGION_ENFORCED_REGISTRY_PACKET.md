# P1175 Safe-Region Enforced Registry Packet

Status: `P1175_EXECUTED_SAFE_REGION_ENFORCED_REGISTRY_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać rekomendowany krok po `P1174`: podpiąć safe-region prefilter jako
opcjonalny twardy warunek w rejestrze kandydatów.

## Professor-level decision

Rozszerzam `P1152` o flagę:

```text
--require-safe-region
```

Semantyka:

1. kandydat poza `P1172 safe_bounds` jest odrzucany przed `P1151`,
2. odrzucenie jest jawnie raportowane (`returncode=9`, `filtered_by_safe_region`),
3. tryb może działać z `--allow-failures` dla raportów mieszanych.

## Demonstration result

Na mieszanym registry (`inside` + `outside`):

- `total=2`, `passed=1`, `failed=1`,
- kandydat outside został odfiltrowany na etapie safe-region.

## Artifacts

- updated runner:
  `p1152_strict_candidate_registry_runner.py`
- sample mixed registry:
  `generated/p1175_mixed_registry.json`
- run summary:
  `generated/p1152_strict_candidate_registry_runner_summary.json`

## Honest boundary

`P1175` wzmacnia rygor selekcji wejścia i nie stanowi claimu closure ani
`QW-2191` discharge.
