# P1213 External CAS Certificate Contract Packet

Status: `P1213_EXECUTED_EXTERNAL_CERTIFICATE_CONTRACT_GATE_NO_FALSE_PASS`
As of: `2026-05-11`

## Goal

Następny uczciwy krok po `P1212`: wymusić minimalny kontrakt certyfikatu
zewnętrznego CAS zanim artefakt zostanie uznany za gotowy do prób W4.

## Professor-level decision

Dodajemy `p1213_external_cas_certificate_contract.py`, który sprawdza, czy
artefakt `P1205` zawiera jawne dane reprodukowalności:

1. `cas_engine` i `cas_engine_version`,
2. `deterministic_replay_seed`,
3. `trace_schema_version`,
4. `trace_hash_sha256` zgodny z formatem,
5. spójność daty i stały status `strict_closure_claim_allowed=false`.

## Honest boundary

`P1213` nie dowodzi W4 i nie tworzy closure claim. To wyłącznie kontrakt
minimalny dla audytu źródła/replay, żeby uniknąć fałszywych przejść bramki.
