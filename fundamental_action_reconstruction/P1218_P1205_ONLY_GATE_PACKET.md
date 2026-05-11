# P1218 P1205-Only Gate Packet

Status: `P1218_EXECUTED_P1205_ONLY_OPERATIONAL_GATE_NO_FALSE_PASS`
As of: `2026-05-11`

## Goal

Uprościć ścieżkę operacyjną: dla wejścia do próby W4 ma wystarczyć `P1205`.

## Professor-level decision

Dodaję w `P1209` tryb:

- `--p1205-only`.

W tym trybie gate opiera się wyłącznie na polach z `P1205`:

1. `sympy_available=true`,
2. `trace_payload` obecny,
3. `trace_hash_sha256` obecny.

## Honest boundary

Tryb `P1205-only` jest tylko uproszczeniem operacyjnym.
Nie oznacza dowodu W4 ani domknięcia teorii (`theory_closure_status=OPEN`).
