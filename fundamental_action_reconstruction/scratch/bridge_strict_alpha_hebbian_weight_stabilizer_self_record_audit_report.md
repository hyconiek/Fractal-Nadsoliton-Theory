# Scratch strict-alpha Hebbian weight stabilizer self-record audit probe

Status: d5 Hebbian weight carries the required {1,11} subgroup conditionally; not a d5 origin theorem.

- D5 unit stabilizer: `[1, 11]`.
- D5 stabilizer equals required subgroup: `True`.
- D5 folded-distance weights: `{'1': '-25/12', '2': '11/12', '3': '-1/12', '4': '-13/12', '5': '23/12', '6': '-25/12'}`.
- Classes with required [1, 11] stabilizer: `['contiguous_step_1_11', 'fifth_step_d5_step_5_7']`.
- Classes with full Aut stabilizer: `['fourth_step_4_8_degenerate', 'nyquist_step_6', 'parity_minus_one_step_2_10', 'third_step_3_9_degenerate']`.
- Contiguous has same stabilizer as d5: `True`.
- Unit 5 obstruction example: `{'source_folded_distance': 1, 'image_folded_distance': 5, 'source_weight': '-25/12', 'image_weight': '23/12'}`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: learned d5 Hebbian weight can carry the subgroup record, but the d5 teacher trace and Hebbian law remain supplied inputs.
- No false pass: no strict d5-origin theorem, no QW-2191 discharge, no ToE closure.
