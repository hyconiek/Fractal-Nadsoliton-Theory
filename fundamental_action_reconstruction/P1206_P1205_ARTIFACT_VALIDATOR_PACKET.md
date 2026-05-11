# P1206 P1205 Artifact Validator Packet

Status: `P1206_EXECUTED_P1205_ARTIFACT_VALIDATOR_NO_FALSE_PASS`
As of: `2026-05-11`

## Goal

Po pobraniu `p1205_w4_sympy_cas_runner_summary.json` dodać formalną bramkę,
która waliduje artefakt przed wejściem do dalszych kroków W4.

## Professor-level decision

Dodaję `p1206_p1205_artifact_validator.py`, który sprawdza:

1. minimalny schemat artefaktu,
2. czy CAS był realnie dostępny,
3. czy istnieją `trace_payload` i `trace_hash_sha256` wymagane do audytu.

## Honest boundary

`P1206` nie wykonuje discharge W4; tylko kontroluje jakość i gotowość
artefaktu wejściowego.
