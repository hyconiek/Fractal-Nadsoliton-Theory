# P1478 — S4.28 QW-2191 SP1 Conservative Policy Update (PL)

Status: `P1478_EXECUTED_QW2191_SP1_POLICY_UPDATE_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Zaktualizować lokalną politykę operacyjną SP1 na podstawie najgorszej granicy
z P1477, z konserwatywnym marginesem bezpieczeństwa.

## Reguła

- `edge_worst = first_fail_shift` z P1477,
- `safe_min_shift = edge_worst + margin`,
- reruny poza `[safe_min_shift, 0]` są blokowane.

## Rygor

- bez legacy bridge,
- bez strict-core closure claim,
- local-only policy artifact.
