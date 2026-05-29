# Scratch strict-alpha Hebbian cycle-metric subgroup selector audit probe

Status: cycle metric supplies the required {1,11} subgroup conditionally; not a strict source theorem.

- Supports scanned: `792` five-active-node states on `Z_12`.
- Distance stabilizers: `{'1': [1, 11], '2': [1, 5, 7, 11], '3': [1, 5, 7, 11], '4': [1, 5, 7, 11], '5': [1, 11], '6': [1, 5, 7, 11]}`.
- Distance `1` stabilizer equals required subgroup: `True`.
- Distance `5` stabilizer equals required subgroup: `True`.
- Full Aut d5 selector invariant: `False`.
- Cycle-metric reduced d5 selector invariant: `True`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: a cycle metric/locality record would supply the needed subgroup, but that record is not derived here.
- No false pass: no strict d5-source theorem, no QW-2191 discharge, no ToE closure.
