# P1483 — S4.33 QW-2191 Strict Sector Split Probe (PL)

Status: `P1483_EXECUTED_QW2191_STRICT_SECTOR_SPLIT_PROBE_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Wykonać pierwszy jawny, liczbowo audytowalny krok konstrukcyjny w torze:

`F(nadsoliton) => L_SM + L_GR`

bez legacy bridge i bez claimu strict-core closure.

## Decyzja profesorska

Budujemy **kandydata splitu sektorowego** z jawnie kontrolowanym członem mieszanym:

- `L_total = w_SM*L_SM + w_GR*L_GR + eps_mix*L_mix`,
- `w_SM + w_GR = 1`, `0 <= eps_mix <= eps_cap`,
- brak ukrytego selektora: `QW-2191` pozostaje otwarte.

## Kryterium uczciwości

Krok uznajemy za uczciwy tylko jeśli:

1. strict-only (`no_legacy_bridge = true`),
2. SP1 local i cross-scenario są dodatnie/przechodzą,
3. `eps_mix` pozostaje małe i jawnie limitowane policy-gate,
4. raport końcowy zawiera jawne `next_step_recommendation` i `layman_explanation`.
