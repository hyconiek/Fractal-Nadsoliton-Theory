# P1210 Controlled W4 Symbolic Attempt Packet

Status: `P1210_EXECUTED_CONTROLLED_W4_SYMBOLIC_ATTEMPT_NO_FALSE_PASS`
As of: `2026-05-11`

## Goal

Po otwarciu `P1209` wykonać kontrolowaną próbę symboliczną W4 z twardą
polityką: discharge tylko przy kompletnym dowodzie (`trace_payload` + hash +
`reduced_zero_ok=true`).

## Professor-level decision

Dodaję `p1210_controlled_w4_symbolic_attempt.py`, który rozdziela:

1. `gate_open`,
2. `attempt_executed`,
3. `w4_discharge_pass`.

## Honest boundary

Nawet jeśli `W4` zostanie oznaczone jako `DISCHARGED`, `strict_closure_claim_allowed`
pozostaje `false` bez dalszych świadectw globalnych.
