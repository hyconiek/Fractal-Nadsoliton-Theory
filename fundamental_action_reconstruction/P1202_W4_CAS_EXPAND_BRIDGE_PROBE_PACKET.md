# P1202 W4 CAS Expand Bridge Probe Packet

Status: `P1202_EXECUTED_W4_CAS_EXPAND_BRIDGE_PROBE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Następny uczciwy krok po `P1201`: sprawdzić pierwszy realny most do backendu
CAS dla operacji `expand`.

## Professor-level decision

Dodaję `p1202_w4_cas_expand_bridge_probe.py`, który:

1. testuje dostępność `sympy`,
2. wykonuje minimalny `expand` na wyrażeniu kontrolnym,
3. raportuje wynik jako `partial_symbolic_execution_ready` bez claimu discharge.

## Honest boundary

`P1202` nie discharge'uje W4; to tylko pierwszy test mostu CAS.
