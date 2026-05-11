# P1196 W4 Symbolic Ledger Template Packet

Status: `P1196_EXECUTED_W4_SYMBOLIC_LEDGER_TEMPLATE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać kolejny uczciwy krok po `P1195`: wyeksportować minimalny template
ledgera symbolicznego dla `W4`, który narzuca format dowodu krok-po-kroku.

## Professor-level decision

Dodaję `p1196_w4_symbolic_ledger_template.py`, który ustala:

1. wymagany nośnik (`psi4**2/2`),
2. sekwencję operacji symbolicznych (expand -> group -> cancel -> reduce),
3. twardy warunek końcowy: `reduced_form == 0`.

## Current outcome

Template został wyeksportowany (`template_exported = true`), ale
`can_discharge_w4_now = false`.

## Honest boundary

`P1196` nie jest dowodem; to formalizacja formatu dowodu.
