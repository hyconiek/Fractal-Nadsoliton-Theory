# P1169 E2E Candidate Integration Packet

Status: `P1169_EXECUTED_E2E_CANDIDATE_INTEGRATION_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać rekomendowany krok po `P1168`: zintegrować zwycięski kandydat
(`sigma=1.2, kappa=0.12`) z pełnym pipeline i registry bez `force_fail_stage`.

## Professor-level decision

Dodaję formalny test E2E:

1. uruchom `P1151` na kandydacie top z `P1166`,
2. uruchom `P1152` na registry z tym kandydatem (bez allow-failures),
3. zapisz wspólny werdykt integracyjny `integrated_pass`.

## Artifacts

- candidate:
  `generated/p1169_candidate_top_from_p1166.json`
- registry:
  `generated/p1169_candidate_registry_top.json`
- integration script:
  `p1169_e2e_candidate_integration.py`
- summary:
  `generated/p1169_e2e_candidate_integration_summary.json`

## Result

`integrated_pass = true` na current repo state.

Czyli kandydat przechodzi pełny workflow metodologiczny jako realny obiekt
roboczy.

## Honest boundary

To nadal walidacja procesu i kandydata proxy. Brak claimu closure i brak
`QW-2191` discharge.


Dodano tryb `--require-shortlist-consistency`, który wymaga w `P1169` aby
`P1152` był uruchomiony z `--rank-by-safe-margin --enforce-shortlist-consistency`
i aby checker `P1177` zwrócił success przed uznaniem integracji za pass.

Dodano tryb `--require-safe-region`, który propaguje filtr bezpiecznego
regionu do `P1152`; w połączeniu z `--require-shortlist-consistency` tworzy to
pełny gate: pipeline-pass + safe-region + shortlist-integrity.

Dodano tryb `--strict-e2e`: wymusza jednocześnie `--require-safe-region` i
`--require-shortlist-consistency`, a dodatkowo wymaga pełnego passu rejestru
(`passed == total`) przed `integrated_pass=true`.

Dodano tryb `--require-out-of-locality-robustness`: `P1169` uruchamia `P1171`
i wymaga wysokiej odporności out-of-locality (`robust_fraction >= 0.99`) przed
`integrated_pass=true`.

Dodano parametr `--robustness-threshold <x>` dla trybu robustness, aby jawnie
ustawiać wymagany próg `robust_fraction` (domyślnie `0.99`) i raportować go w
artefakcie `P1169`.
