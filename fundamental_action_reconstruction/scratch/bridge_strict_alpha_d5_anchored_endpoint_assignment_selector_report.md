# Scratch strict-alpha d5 anchored endpoint assignment selector probe

Status: anchored endpoint assignment selector; no strict source theorem.

- Reference ordered support `[0, 7, 2, 9, 4]` from source `0`, orientation `-1`.
- Before endpoint bias: path minimizer count `2`.
- Forward endpoint selector chooses `[2, 2, 2, 1, 1]` uniquely: `True`.
- Orbit scan: `24` ordered supports; all forward assignments unique: `True`.
- Target replay: `q^5=256/243`, eta residual `0.000e+00`.
- Honest read: source+orientation+endpoint bias closes the finite assignment table, but those are extra selector premises.
- No false pass: no strict source/orientation/action theorem, no QW-2191 discharge, no ToE closure.
