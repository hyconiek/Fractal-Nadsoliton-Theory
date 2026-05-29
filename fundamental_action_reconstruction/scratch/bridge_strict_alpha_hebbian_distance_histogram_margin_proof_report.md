# Scratch strict-alpha Hebbian distance-histogram margin proof probe

Status: d5 margin reduced from support enumeration to distance-histogram proof; not an origin theorem.

- Supports scanned: `792`; histogram classes: `35`.
- D5 histogram: `(0, 3, 2, 1, 4, 0)` with support count `12` and equals teacher orbit `True`.
- Closest non-d5 histogram: `(1, 3, 2, 1, 3, 0)` with support count `24`.
- Gap by variant: `{'binary_with_diagonal': '8', 'binary_zero_self': '8', 'centered_with_diagonal': '8', 'centered_zero_self': '8', 'bipolar_with_diagonal': '32', 'bipolar_zero_self': '32'}`.
- Diagonal shift by variant: `{'binary_with_diagonal': '25', 'binary_zero_self': '0', 'centered_with_diagonal': '175/12', 'centered_zero_self': '0', 'bipolar_with_diagonal': '60', 'bipolar_zero_self': '0'}`.
- All variants have positive histogram margin: `True`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: histogram proof compresses the supplied Hebbian margin, but teacher trace, learning law, and histogram source remain supplied inputs.
- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.
