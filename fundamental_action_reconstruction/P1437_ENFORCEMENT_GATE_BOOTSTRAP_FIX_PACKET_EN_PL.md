# P1437 — Enforcement Gate Bootstrap Fix

Status: `P1437_GATE_BOOTSTRAP_EDGE_FIXED`
As of: `2026-05-13`

## Professorial decision

Fix the self-reference bootstrap edge in P1435 gate so the checker does not fail on its own output file.

## Fix applied

1. Exclude `p1435_pass_scope_enforcement_gate_summary.json` from audited target set.
2. Ensure P1435 summary itself writes mandatory semantics fields.

## Result

Gate can now evaluate strict-route summaries without self-triggered false violation.

## PL (skrót)

Naprawiono problem "kontroler sprawdza samego siebie".
