# P1496 — S4.46 QW-2191 Cross-Provider Replication Packet (PL)

Status: `P1496_EXECUTED_QW2191_CROSS_PROVIDER_REPLICATION_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Wykonać niezależną replikację kandydata twierdzenia (`P1495`) na dwóch
odseparowanych klasach dostawcy danych (provider class A/B) bez legacy bridge.

## Decyzja profesorska

Dzielimy bezpieczny zakres na dwa niezależne zbiory:

- Provider A: punkty o indeksie parzystym,
- Provider B: punkty o indeksie nieparzystym.

W obu zbiorach musi zajść:

1. `G1 < G0`,
2. `|Delta_SB| <= margin`,
3. zgodna orientacja selektora.

Jeśli oba przechodzą, mamy lokalną replikację cross-provider.
