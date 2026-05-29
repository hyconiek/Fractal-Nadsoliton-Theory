# Scratch strict-alpha phase-mode aliasing obstruction probe

Status: phase-mode aliasing classification under phase-origin premise; no strict selector discharge.

- Source-complete modes: `[1, 5, 7, 11]`; aliasing modes: `[0, 2, 3, 4, 6, 8, 9, 10]`.
- Expected phase class counts by mode: `{'0': 1, '1': 12, '2': 6, '3': 4, '4': 3, '5': 12, '6': 2, '7': 12, '8': 3, '9': 4, '10': 6, '11': 12}`.
- All observed counts match gcd formula: `True`.
- Target replay: `q^5=256/243`, eta residual `0.000e+00`.
- Honest read: a single phase mode recovers source iff gcd(k,12)=1; non-coprime modes alias sources.
- No false pass: no strict phase-origin theorem, no QW-2191 discharge, no ToE closure.
