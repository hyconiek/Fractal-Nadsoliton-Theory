# Scratch strict-alpha existence bit vs unit-orientation selector audit

Status: existence bit gates nonempty orbit but does not select d5 branch.

- Supports scanned: `792`; histogram classes: `35`.
- Existence gate functions audited: `64`.
- Canonical existence gates E=0->empty, E=1->{A1,A5}: `1`.
- Full-Aut invariant singleton-d5 gates from E=1: `0`.
- Direct correction: existence/nonexistence is an Aut-invariant scalar bit; unit orientation is an Aut-breaking bit.
- Honest equation: `E=0 -> empty`, `E=1 -> {A1,A5}`.
- Forbidden shortcut: `E=1 -> A5/d5` silently imports the unit-orientation bit.
- Target replay kept conditional: `q^5=256/243`, eta `9/5`.
- No false pass: no existence-bit source theorem, no QW-2191 discharge, no ToE closure.
