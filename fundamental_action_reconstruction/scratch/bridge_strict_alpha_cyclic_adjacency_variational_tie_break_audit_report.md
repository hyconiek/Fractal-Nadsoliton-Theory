# Scratch strict-alpha cyclic-adjacency variational tie-break audit probe

Status: nearest-neighbour adjacency alone is not a global d5 selector; d1/d5 tie-break is conditional.

- Supports scanned: `792`; histogram classes: `35`.
- A1 histogram d1..d6: `[4, 3, 2, 1, 0, 0]`.
- A5 histogram d1..d6: `[0, 3, 2, 1, 4, 0]`.
- Minimum d1 support count: `36`.
- Minimum d1 dihedral-orbit count: `3`.
- d1 alone uniquely selects A5: `False`.
- Lexicographic min d1 / max d5 selects A5: `True`.
- Lexicographic winner support count: `12`; orbit count: `1`.
- Target replay kept conditional: `q^5=256/243`, eta `9/5`.
- No false pass: no strict variational-source theorem, no QW-2191 discharge, no ToE closure.
