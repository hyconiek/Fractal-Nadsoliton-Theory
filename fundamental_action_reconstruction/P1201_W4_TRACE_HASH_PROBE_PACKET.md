# P1201 W4 Trace-Hash Probe Packet

Status: `P1201_EXECUTED_W4_TRACE_HASH_PROBE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Następny uczciwy krok po `P1200`: uruchomić audytowalność trace'u przez
stabilny hash, nawet jeśli backend symboliczny jeszcze nie wykonuje algebra.

## Professor-level decision

Dodaję `p1201_w4_trace_hash_probe.py`, który eksportuje:

1. standaryzowany payload trace'u,
2. hash SHA-256 do audytu reprodukowalności,
3. rozdzielenie: `auditability_ready=true` vs `symbolic_execution_ready=false`.

## Current outcome

Audytowalność techniczna jest gotowa, ale wykonanie symboliczne nadal nie.

## Honest boundary

`P1201` nie jest discharge W4; to krok infrastruktury audytu.
