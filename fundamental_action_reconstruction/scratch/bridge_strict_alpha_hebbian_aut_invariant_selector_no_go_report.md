# Scratch strict-alpha Hebbian Aut-invariant selector no-go probe

Status: finite Aut-invariant selector no-go for the d5 unit pair.

- Supports scanned: `792` five-active-node states on `Z_12`.
- Aut action table: `{'1': {'contiguous_step_1_11': 'contiguous_step_1_11', 'parity_minus_one_step_2_10': 'parity_minus_one_step_2_10', 'fifth_step_d5_step_5_7': 'fifth_step_d5_step_5_7'}, '5': {'contiguous_step_1_11': 'fifth_step_d5_step_5_7', 'parity_minus_one_step_2_10': 'parity_minus_one_step_2_10', 'fifth_step_d5_step_5_7': 'contiguous_step_1_11'}, '7': {'contiguous_step_1_11': 'fifth_step_d5_step_5_7', 'parity_minus_one_step_2_10': 'parity_minus_one_step_2_10', 'fifth_step_d5_step_5_7': 'contiguous_step_1_11'}, '11': {'contiguous_step_1_11': 'contiguous_step_1_11', 'parity_minus_one_step_2_10': 'parity_minus_one_step_2_10', 'fifth_step_d5_step_5_7': 'fifth_step_d5_step_5_7'}}`.
- Orbit partition: `[['contiguous_step_1_11', 'fifth_step_d5_step_5_7'], ['parity_minus_one_step_2_10']]`.
- Invariant subsets: `[[], ['parity_minus_one_step_2_10'], ['contiguous_step_1_11', 'fifth_step_d5_step_5_7'], ['contiguous_step_1_11', 'parity_minus_one_step_2_10', 'fifth_step_d5_step_5_7']]`.
- d5 singleton Aut-invariant: `False`.
- Unit pair Aut-invariant: `True`.
- Previous selector classification: `symmetry-breaking selector premise`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: strict Aut-invariant logic can keep {contiguous,d5} together but cannot choose d5 alone.
- No false pass: no strict d5-source theorem, no QW-2191 discharge, no ToE closure.
