# Scratch strict-alpha Hebbian symbolic mode-competition no-go probe

Status: unconstrained highest resonance selects Nyquist k=6, not d5.

- Supports scanned: `792`; histogram classes: `35`.
- k=5 top: `{'distance_histogram_d1_to_d6': [0, 3, 2, 1, 4, 0], 'support_count': 12, 'power_exact': {'rational_part': '7', 'sqrt3_coefficient': '4', 'expression': '7 + 4*sqrt3'}}`.
- k=6 top: `{'distance_histogram_d1_to_d6': [0, 4, 0, 4, 0, 2], 'support_count': 12, 'power_exact': {'rational_part': '25', 'sqrt3_coefficient': '0', 'expression': '25'}}`.
- Global winner mode: `6` with power `25`.
- Fifth-mode d5 power: `7 + 4*sqrt3`; global-minus-fifth `18 - 4*sqrt3`.
- Unconstrained highest resonance selects d5: `False`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: k=5 selects d5 only conditionally; unconstrained mode competition selects k=6/Nyquist.
- No false pass: no strict fifth-mode source theorem, no QW-2191 discharge, no ToE closure.
