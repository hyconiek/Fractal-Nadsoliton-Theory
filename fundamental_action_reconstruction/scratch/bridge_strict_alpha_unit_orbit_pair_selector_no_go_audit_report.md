# Scratch strict-alpha unit-orbit pair selector no-go audit probe

Status: full Aut forbids singleton d5-pair selection; chi_11-kernel reduction allows it conditionally.

- Branch modes audited: `[1, 5, 7, 11]`.
- Full-Aut invariant singleton d5 selectors: `0`.
- chi_11-kernel invariant singleton d5 selectors: `1`.
- Full-Aut invariant selectors: `[{'selector': 'empty_selector', 'selected_pairs': [], 'selects_singleton_d5_pair': False, 'invariant_under': 'full_Aut_Z12', 'is_invariant': True}, {'selector': 'both_pairs_unoriented_orbit', 'selected_pairs': ['contiguous_pair_A1_A11', 'd5_pair_A5_A7'], 'selects_singleton_d5_pair': False, 'invariant_under': 'full_Aut_Z12', 'is_invariant': True}]`.
- chi_11-kernel invariant selectors: `[{'selector': 'empty_selector', 'selected_pairs': [], 'selects_singleton_d5_pair': False, 'invariant_under': 'chi_11_kernel_units_{1,11}', 'is_invariant': True}, {'selector': 'singleton_contiguous_pair', 'selected_pairs': ['contiguous_pair_A1_A11'], 'selects_singleton_d5_pair': False, 'invariant_under': 'chi_11_kernel_units_{1,11}', 'is_invariant': True}, {'selector': 'singleton_d5_pair', 'selected_pairs': ['d5_pair_A5_A7'], 'selects_singleton_d5_pair': True, 'invariant_under': 'chi_11_kernel_units_{1,11}', 'is_invariant': True}, {'selector': 'both_pairs_unoriented_orbit', 'selected_pairs': ['contiguous_pair_A1_A11', 'd5_pair_A5_A7'], 'selects_singleton_d5_pair': False, 'invariant_under': 'chi_11_kernel_units_{1,11}', 'is_invariant': True}]`.
- Target replay kept conditional: `q^5=256/243`, eta `9/5`.
- No false pass: no strict symmetry-reduction source theorem, no QW-2191 discharge, no ToE closure.
