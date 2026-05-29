# Scratch strict-alpha Hebbian rule-variant robustness audit probe

Status: d5 Hebbian stabilizer and energy maxima are robust across tested readouts; not an origin theorem.

- Tested variants: `6`.
- Passing variants: `['binary_with_diagonal', 'binary_zero_self', 'centered_with_diagonal', 'centered_zero_self', 'bipolar_with_diagonal', 'bipolar_zero_self']`.
- Failing variants: `[]`.
- All tested variants have required stabilizer: `True`.
- All tested variants have unique d5 global maxima: `True`.
- Maximum energy by variant: `{'binary_with_diagonal': '85', 'binary_zero_self': '60', 'centered_with_diagonal': '395/12', 'centered_zero_self': '55/3', 'bipolar_with_diagonal': '140', 'bipolar_zero_self': '80'}`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: robust inside tested Hebbian-family readouts, but teacher trace and learning law remain supplied inputs.
- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.
