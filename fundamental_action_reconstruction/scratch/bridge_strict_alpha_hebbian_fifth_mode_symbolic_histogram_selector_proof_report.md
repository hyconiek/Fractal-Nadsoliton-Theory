# Scratch strict-alpha Hebbian fifth-mode symbolic histogram selector proof probe

Status: fifth-mode power selects d5 symbolically; not an origin theorem.

- Supports scanned: `792`; histogram classes: `35`.
- P5 formula: `P_5(h)=5+2*sum_d h_d*cos(5*pi*d/6)` in `Q(sqrt(3))`.
- Top histogram: `{'distance_histogram_d1_to_d6': [0, 3, 2, 1, 4, 0], 'support_count': 12, 'p5_power_exact': {'rational_part': '7', 'sqrt3_coefficient': '4', 'expression': '7 + 4*sqrt3'}}`.
- Second histogram: `{'distance_histogram_d1_to_d6': [1, 3, 2, 1, 3, 0], 'support_count': 24, 'p5_power_exact': {'rational_part': '7', 'sqrt3_coefficient': '2', 'expression': '7 + 2*sqrt3'}}`.
- Symbolic gap: `2*sqrt3`; positive `True`.
- Fifth-mode channel stabilizer: `[1, 11]`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: exact fifth-mode readout selects d5, but k=5 remains a supplied selector channel.
- No false pass: no strict fifth-mode source theorem, no QW-2191 discharge, no ToE closure.
