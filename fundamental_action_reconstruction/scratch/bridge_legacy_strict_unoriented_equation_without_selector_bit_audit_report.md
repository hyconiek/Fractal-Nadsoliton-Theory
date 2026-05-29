# Scratch legacy->strict unoriented bridge equation without selector bit audit

Status: quotient bridge equation possible; singleton d5 bridge blocked without unit bit.

- Supports scanned: `792`; histogram classes: `35`.
- Full-Aut invariant bridge outputs: `['empty_relation', 'unoriented_orbit_A1_A5']`.
- Singleton d5 full-Aut invariant: `False`.
- Unoriented orbit output available: `True`.
- Direct answer: without the one-bit selector, only `Bridge_0(K_legacy_ont) -> {A1,A5}` is honest.
- Forbidden shortcut: `K_legacy_ont -> K_strict_gate(d5)` would silently import the selector bit.
- Target replay kept as conditional data: `q^5=256/243`, eta `9/5`.
- No false pass: no K_legacy_ont == K_strict_gate theorem, no QW-2191 discharge, no ToE closure.
