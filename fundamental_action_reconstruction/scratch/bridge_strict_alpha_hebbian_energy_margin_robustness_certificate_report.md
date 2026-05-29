# Scratch strict-alpha Hebbian energy-margin robustness certificate probe

Status: positive d5 energy margins across tested Hebbian readouts; not an origin theorem.

- Tested variants: `6`.
- Gap by variant: `{'binary_with_diagonal': '8', 'binary_zero_self': '8', 'centered_with_diagonal': '8', 'centered_zero_self': '8', 'bipolar_with_diagonal': '32', 'bipolar_zero_self': '32'}`.
- Minimum gap: `8`.
- Best non-d5 competitor count by variant: `{'binary_with_diagonal': 24, 'binary_zero_self': 24, 'centered_with_diagonal': 24, 'centered_zero_self': 24, 'bipolar_with_diagonal': 24, 'bipolar_zero_self': 24}`.
- Support-energy perturbation safe radius by variant: `{'binary_with_diagonal': '4', 'binary_zero_self': '4', 'centered_with_diagonal': '4', 'centered_zero_self': '4', 'bipolar_with_diagonal': '16', 'bipolar_zero_self': '16'}`.
- Entrywise-weight perturbation sufficient bound by variant: `{'binary_with_diagonal': '4/25', 'binary_zero_self': '4/25', 'centered_with_diagonal': '4/25', 'centered_zero_self': '4/25', 'bipolar_with_diagonal': '16/25', 'bipolar_zero_self': '16/25'}`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: exact positive margins protect the supplied Hebbian readouts, but the teacher trace, learning law, and perturbation norm remain supplied inputs.
- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.
