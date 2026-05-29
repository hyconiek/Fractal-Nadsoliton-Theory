# Scratch strict-alpha Hebbian combinatorial score margin proof probe

Status: d5 margin reduced to an integer histogram-score inequality; not an origin theorem.

- Supports scanned: `792`; histogram classes: `35`.
- Score formula: `C(h)=3*h_2 + 2*h_3 + h_4 + 4*h_5`; coefficients `(0, 3, 2, 1, 4, 0)`.
- Top histogram row: `{'distance_histogram_d1_to_d6': [0, 3, 2, 1, 4, 0], 'support_count': 12, 'combinatorial_score': 30}`.
- Second histogram row: `{'distance_histogram_d1_to_d6': [1, 3, 2, 1, 3, 0], 'support_count': 24, 'combinatorial_score': 26}`.
- Score gap: `4`; variant gap replay: `{'binary_with_diagonal': 8, 'binary_zero_self': 8, 'centered_with_diagonal': 8, 'centered_zero_self': 8, 'bipolar_with_diagonal': 32, 'bipolar_zero_self': 32}`.
- D5 histogram supports equal teacher orbit: `True`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: integer score proves the supplied Hebbian margin, but score, teacher trace, and learning law remain supplied inputs.
- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.
